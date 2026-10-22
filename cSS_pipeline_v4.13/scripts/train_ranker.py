#!/usr/bin/env python3
"""Fit the candidate ranking model on expert labels (v5.0, Layer 3) -> review/ranker_model.json.
usage: train_ranker.py [--split dev] [--source consensus|reader:NAME] [--l2 1.0] [--min-patients 3] [--boot 500]
                       [--out review/ranker_model.json]
What it does:
  - training rows = shown candidates AND excluded objects (from expert sheets made with --include-excluded) that
    carry an expert label in labels/consensus.csv (or one reader's latest session in candidate_labels.csv),
    matched by cand_uid fingerprint, from studies of the given split ONLY (labels/study_meta.csv). The default
    split is dev, and --split test is refused: test studies are for reporting with evaluate.py.
  - target: 1 = called cSS, 0 = called anything else; 'uncertain' rows are left out.
  - features: FEATURES below, read from review/ID_candidates|excluded.csv plus review/ID_surface_features.csv
    (recon-all route) when present. A feature missing from every training row is dropped with a warning;
    missing values are filled with the training median. Deprecated columns (vein_tree_mm, mirror_dark_frac,
    branch_per10mm) and the scores (score, score_v3, score_v4) are NOT features.
  - model: L2-regularised logistic regression (scipy.optimize.minimize, L-BFGS), class-balanced weights,
    features median-imputed and standardised. Leave-one-patient-out (LOPO) cross-validation: imputation and
    standardisation are refitted inside every fold. Reported: pooled LOPO AUC of the model AND of score_v4 on
    the same rows (the baseline it must beat), each with a patient-bootstrap 95 % interval, and per-patient
    numbers. The final model is then fitted on all training rows.
Needs >= --min-patients labelled patients with >= 1 positive and >= 1 negative row overall; stops otherwise.
With fewer than ~10 patients every number here is noisy: read the intervals, not the point estimates.
rank_candidates.py applies the model. Research software, not for clinical use."""
import sys, os, json, argparse, numpy as np, pandas as pd
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, cand_uid, label_class, FP

FEATURES = [
    # darkness and pial position
    "darkness_z", "edge_contrast", "pial_dist_mm", "surface_contact", "cortex_frac", "cortex_dist_mm", "bank_frac",
    # shape
    "volume_mm3", "n_slices", "elongation", "extent_mm", "tube_ratio", "surface_alignment", "tram_frac",
    "parenchyma_frac", "n_pieces",
    # vessel continuation (bounded)
    "vessel_wm_mm", "vessel_ext_mm", "vessel_run_mm", "axis_normal", "local_cont_max_mm", "local_cont_min_mm",
    "terminal_ends", "local_branch_per10mm", "net_radius_med_mm",
    # thickness
    "cand_radius_med_mm", "cand_radius_p95_mm", "cand_radius_cv",
    # regional / contralateral reference
    "region_z", "region_dark_frac", "mirror_region_z", "cand_mirror_z",
    # FLAIR and mimic evidence
    "flair_csf_z", "flair_ctx_z", "mimic_score", "phase_z", "mb_cluster_n",
    # pial surface (recon-all route, review/ID_surface_features.csv)
    "surf_normal_dist_mm", "surf_width_mm", "surf_inradius_mm", "surf_aspect", "geo_length_mm", "surf_band_continuity",
    "surf_dark_cont_mm", "surf_exit_paths", "surf_sulc", "surf_curv", "surf_thickness_mm"]
DEPRECATED = ["vein_tree_mm", "mirror_dark_frac", "branch_per10mm"]       # never features (saturated / uninformative)


def load_labels(L, source):
    if source == "consensus":
        lab = pd.read_csv(f"{L}/consensus.csv"); lab = lab[lab.method != "pending"]
        return lab.assign(cls=lab.consensus_class)[["study_id", "cand_uid", "cls"]]
    if source.startswith("reader:"):
        rd = source.split(":", 1)[1]; lab = pd.read_csv(f"{L}/candidate_labels.csv"); lab = lab[lab.reader_id == rd]
        lab = lab[lab.session == lab.groupby("study_id").session.transform("max")]
        return lab.assign(cls=lab.label.map(label_class))[["study_id", "cand_uid", "cls"]]
    sys.exit("--source consensus or reader:NAME")


def study_table(B, s):
    """candidates + excluded of one study (with surface features merged on the fingerprint when present)"""
    parts = []
    for lst in ("candidates", "excluded"):
        p = f"{B}/review/{s}_{lst}.csv"
        if os.path.exists(p):
            d = pd.read_csv(p)
            if len(d): parts.append(d.assign(_list=lst))
    if not parts: return None
    d = pd.concat(parts, ignore_index=True)
    sp = f"{B}/review/{s}_surface_features.csv"
    if os.path.exists(sp):
        sf = pd.read_csv(sp).drop(columns=["list", "obj_id"], errors="ignore").drop_duplicates(FP)
        d = d.merge(sf, on=FP, how="left", suffixes=("", "_surf"))
    d["cand_uid"] = [cand_uid(s, r) for _, r in d.iterrows()]
    return d


def prep(X, med=None, mu=None, sd=None):
    """median imputation + standardisation; parameters from X unless given"""
    X = np.asarray(X, float)
    if med is None: med = np.nanmedian(X, axis=0); med = np.where(np.isfinite(med), med, 0.0)
    X = np.where(np.isfinite(X), X, med)
    if mu is None: mu = X.mean(0); sd = X.std(0); sd = np.where(sd > 1e-9, sd, 1.0)
    return (X - mu) / sd, med, mu, sd


def fit_lr(Z, y, l2):
    """class-balanced L2 logistic regression -> (w, b)"""
    w_pos = 0.5 / max(y.mean(), 1e-9); w_neg = 0.5 / max(1 - y.mean(), 1e-9)
    cw = np.where(y == 1, w_pos, w_neg)
    def f(p):
        w, b = p[:-1], p[-1]; t = Z @ w + b
        ll = cw * (np.logaddexp(0, t) - y * t)
        g = cw * (1 / (1 + np.exp(-t)) - y)
        return ll.sum() / len(y) + 0.5 * l2 * (w @ w) / len(y), \
            np.r_[Z.T @ g / len(y) + l2 * w / len(y), g.sum() / len(y)]
    r = minimize(f, np.zeros(Z.shape[1] + 1), jac=True, method="L-BFGS-B")
    return r.x[:-1], float(r.x[-1])


def auc(y, s):
    y = np.asarray(y); s = np.asarray(s, float); pos, neg = s[y == 1], s[y == 0]
    if not len(pos) or not len(neg): return np.nan
    return float(((pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()) / (len(pos) * len(neg)))


def boot_auc(df, col, n, rng):
    pats = df.study_id.unique(); vals = []
    for _ in range(n):
        pick = rng.choice(pats, len(pats), replace=True)
        d = pd.concat([df[df.study_id == p] for p in pick]); vals.append(auc(d.y, d[col]))
    v = np.array([x for x in vals if np.isfinite(x)])
    return [round(float(np.percentile(v, 2.5)), 3), round(float(np.percentile(v, 97.5)), 3)] if len(v) else [None, None]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--split", default="dev", choices=["dev", "all"])   # test is refused: reporting only
    ap.add_argument("--source", default="consensus"); ap.add_argument("--l2", type=float, default=1.0)
    ap.add_argument("--min-patients", type=int, default=3); ap.add_argument("--boot", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--out", default="")
    a = ap.parse_args(); B = _base(); L = f"{B}/labels"
    out = a.out or f"{B}/review/ranker_model.json"
    meta = pd.read_csv(f"{L}/study_meta.csv", dtype=str)
    studies = meta[meta.split == "dev"] if a.split == "dev" else meta[meta.split != "test"]
    if (meta.split == "test").any() and a.split == "all":
        print(f"note: {int((meta.split == 'test').sum())} test-split studies are left out (never trained on)")
    lab = load_labels(L, a.source)
    rows = []
    for s in studies.study_id:
        d = study_table(B, s)
        ls = lab[(lab.study_id == s) & lab.cls.isin([c for c in lab.cls.unique() if c != "uncertain"])]
        if d is None or not len(ls): continue
        m = d.merge(ls[["cand_uid", "cls"]], on="cand_uid", how="inner")
        if len(m) < len(ls): print(f"warning: {s}: {len(ls) - len(m)} labels no longer match the detector output")
        if len(m): rows.append(m.assign(study_id=s, y=(m.cls == "cSS").astype(int)))
    if not rows: sys.exit("no labelled candidates in the requested split - label some studies first (expert_sheet.py)")
    T = pd.concat(rows, ignore_index=True)
    n_pat = T.study_id.nunique()
    print(f"train_ranker: {len(T)} labelled objects ({int(T.y.sum())} cSS, {int((1 - T.y).sum())} other) "
          f"from {n_pat} patients, split {a.split}, source {a.source}")
    if n_pat < a.min_patients: sys.exit(f"need >= {a.min_patients} labelled patients, have {n_pat} - not fitted")
    if T.y.sum() == 0 or T.y.sum() == len(T): sys.exit("need both cSS and non-cSS labels - not fitted")
    feats = [f for f in FEATURES if f in T.columns and T[f].notna().any()]
    gone = [f for f in FEATURES if f not in feats]
    if gone: print(f"warning: {len(gone)} features not available in any training row, dropped: {', '.join(gone)}")
    X = T[feats].apply(pd.to_numeric, errors="coerce").values.astype(float); y = T.y.values
    # leave-one-patient-out: imputation + scaling refitted inside every fold
    T["score_lopo"] = np.nan
    for p in T.study_id.unique():
        tr, te = (T.study_id != p).values, (T.study_id == p).values
        if y[tr].min() == y[tr].max(): continue
        Ztr, med, mu, sd = prep(X[tr]); w, b = fit_lr(Ztr, y[tr], a.l2)
        T.loc[te, "score_lopo"] = 1 / (1 + np.exp(-(prep(X[te], med, mu, sd)[0] @ w + b)))
    base_col = "score_v4" if "score_v4" in T.columns else "score"
    rng = np.random.default_rng(a.seed); ok = T.score_lopo.notna()
    auc_m, auc_b = auc(T.y[ok], T.score_lopo[ok]), auc(T.y[ok], T[base_col][ok])
    ci_m = boot_auc(T[ok], "score_lopo", a.boot, rng); ci_b = boot_auc(T[ok], base_col, a.boot, rng)
    per = [dict(study_id=p, n=int((T.study_id == p).sum()), n_css=int(T.y[T.study_id == p].sum()),
                auc_model=auc(T.y[T.study_id == p], T.score_lopo[T.study_id == p]),
                auc_score_v4=auc(T.y[T.study_id == p], T[base_col][T.study_id == p])) for p in T.study_id.unique()]
    Z, med, mu, sd = prep(X); w, b = fit_lr(Z, y, a.l2)
    model = dict(model="logistic_l2_balanced", version="v5.0", fitted_at=pd.Timestamp.now().isoformat(timespec="seconds"),
                 features=feats, weights=[round(float(v), 6) for v in w], bias=round(b, 6),
                 impute_median=[float(v) for v in med], mean=[float(v) for v in mu], sd=[float(v) for v in sd],
                 l2=a.l2, split=a.split, source=a.source, training_studies=sorted(T.study_id.unique().tolist()),
                 n_objects=int(len(T)), n_css=int(T.y.sum()), n_patients=int(n_pat),
                 lopo_auc=None if not np.isfinite(auc_m) else round(auc_m, 3), lopo_auc_ci95=ci_m,
                 baseline=base_col, baseline_auc=None if not np.isfinite(auc_b) else round(auc_b, 3), baseline_auc_ci95=ci_b,
                 per_patient=[{k: (round(v, 3) if isinstance(v, float) and np.isfinite(v) else (None if isinstance(v, float) else v))
                               for k, v in r.items()} for r in per],
                 dropped_features=gone,
                 note="research software; pooled LOPO AUC over few patients is noisy - read the intervals")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(model, open(out, "w"), indent=2)
    top = sorted(zip(feats, w), key=lambda t: -abs(t[1]))[:8]
    print(f"LOPO AUC model {model['lopo_auc']} {ci_m} vs {base_col} {model['baseline_auc']} {ci_b} (patient bootstrap 95 %)")
    print("largest standardised weights: " + ", ".join(f"{f} {v:+.2f}" for f, v in top))
    print(f"-> {out}   (apply: rank_candidates.py SUBJECT; never trained on test-split studies)")
