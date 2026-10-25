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
  - model: L2-regularised logistic regression (css_model.py, L-BFGS), class-balanced weights, features
    median-imputed and standardised. Leave-one-patient-out (LOPO) cross-validation: imputation and
    standardisation are refitted inside every fold; v5.2: the L2 strength is chosen by an inner LOPO inside each
    training fold (--no-nested: fixed --l2), folds whose training part has one class are skipped and logged.
    Reported on the same folds: the model, score_v4 (the baseline it must beat) and one-family models (darkness
    only, vessel only): AUC with a patient-bootstrap 95 % interval, PR-AUC, sensitivity at 1 and 3 FP per scan,
    Brier (fitted models only); per-patient AUCs. The final model is fitted on all training rows; the JSON
    carries a model hash, the build and research_only; review/ranker_model_card.md summarises it.
Needs >= --min-patients labelled patients with >= 1 positive and >= 1 negative row overall; stops otherwise.
With fewer than ~10 patients every number here is noisy: read the intervals, not the point estimates.
rank_candidates.py applies the model. Research software, not for clinical use."""
import sys, os, json, argparse, numpy as np, pandas as pd
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, cand_uid, label_class, FP, use_derived, sha256, build_id
import css_model as CM

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
    "surf_dark_cont_mm", "surf_exit_paths", "surf_sulc", "surf_curv", "surf_thickness_mm",
    # v5.1.1 missing-input indicators: a missing T1 / FLAIR / phase / surface is information, not a normal value
    "missing_t1_anatomy", "missing_flair_features", "missing_phase_features", "missing_surface_features"]
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
    # v5.1: surface features must postdate this detection (candidates.nii.gz is written only by detect_css.py)
    if use_derived(sp, [f"{B}/review/{s}_candidates.nii.gz"], f"{s} surface features"):
        sf = pd.read_csv(sp).drop(columns=["list", "obj_id"], errors="ignore").drop_duplicates(FP)
        d = d.merge(sf, on=FP, how="left", suffixes=("", "_surf"))
    d["cand_uid"] = [cand_uid(s, r) for _, r in d.iterrows()]
    return add_missing(d)


def add_missing(d):
    """v5.1.1: missing-input indicators (added BEFORE median imputation). detect_css.py writes the T1 / FLAIR / phase
    ones; older outputs and the surface family are derived here from the feature values themselves."""
    d = d.copy()
    surf = [c for c in d.columns if c.startswith("surf_") or c.startswith("geo_")]
    d["missing_surface_features"] = (d[surf].isna().all(axis=1).astype(int) if surf else 1)
    if "missing_phase_features" not in d.columns:
        d["missing_phase_features"] = (d["phase_z"].isna().astype(int) if "phase_z" in d.columns else 1)
    if "missing_flair_features" not in d.columns:
        fl = [c for c in ("flair_csf_z", "flair_ctx_z") if c in d.columns]
        d["missing_flair_features"] = (d[fl].isna().all(axis=1).astype(int) if fl else 1)
    return d


def write_card(M, path):
    """v5.2 model card: what was fitted, on what, how it was evaluated, and what it must not be used for"""
    ev = M["evaluation"]; L = [f"# Candidate ranking model card ({M['version']})", "",
         "Research software. Ranks detector objects for EXPERT REVIEW; not a diagnosis, not a calibrated probability, "
         "not clinically validated. Every report needs expert review.", "",
         f"- fitted: {M['fitted_at']}, build {M['build']}, model hash `{M['model_hash']}`",
         f"- training: {M['n_objects']} labelled objects ({M['n_css']} cSS) from {M['n_patients']} patients, split "
         f"`{M['split']}` (test split never used), labels `{M['source']}` (sha256 {str(M['label_sha256'])[:12]})",
         f"- model: class-balanced L2 logistic regression, {len(M['features'])} features, L2 {M['l2']} ({M['l2_selection']})",
         f"- dropped features (absent in all training rows): {', '.join(M['dropped_features']) or 'none'}",
         "- missing T1 / FLAIR / phase / surface inputs enter as indicator features before median imputation", "",
         "## Leave-one-patient-out evaluation (same folds for every row)", "",
         "| scorer | AUC | 95 % CI (patient bootstrap) | PR-AUC | sens @ 1 FP/scan | sens @ 3 FP/scan | Brier |",
         "|---|---|---|---|---|---|---|"]
    for k, e in ev.items():
        L.append(f"| {k} | {e.get('auc')} | {e.get('auc_ci95')} | {e.get('pr_auc')} | {e.get('sens_at_1fp_per_scan')} | "
                 f"{e.get('sens_at_3fp_per_scan')} | {e.get('brier') if e.get('brier') is not None else '-'} |")
    sk = [f["study_id"] for f in M["lopo_fold_status"] if f["status"] == "skipped"]
    L += ["", f"Folds skipped (training part had one class): {', '.join(sk) or 'none'}.",
          "Brier is given for the fitted models only; score_v4 is a rule score, not a probability.", "",
          "## Limits", "", "- few patients: read the intervals, not the point estimates",
          "- the model must beat score_v4 on held-out TEST studies (evaluate.py) before it changes any reading order by default",
          "- labels from one reader carry that reader's bias; consensus labels are preferred",
          "- not applied to its own training studies (rank_candidates.py refuses them)"]
    open(path, "w").write("\n".join(L) + "\n"); return path


prep, fit_lr, auc = CM.prep, CM.fit_lr, CM.auc          # v5.2: shared with classify_surface_patches.py


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--split", default="dev", choices=["dev", "all"])   # test is refused: reporting only
    ap.add_argument("--source", default="consensus"); ap.add_argument("--l2", type=float, default=1.0)
    ap.add_argument("--min-patients", type=int, default=3); ap.add_argument("--boot", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--out", default="")
    ap.add_argument("--no-nested", action="store_true", help="fixed --l2 instead of nested LOPO selection")
    a = ap.parse_args(); B = _base(); L = f"{B}/labels"
    out = a.out or f"{B}/review/ranker_model.json"
    meta = pd.read_csv(f"{L}/study_meta.csv", dtype=str)
    studies = meta[meta.split == "dev"] if a.split == "dev" else meta[meta.split != "test"]
    if (meta.split == "test").any() and a.split == "all":
        print(f"note: {int((meta.split == 'test').sum())} test-split studies are left out (never trained on)")
    lab = load_labels(L, a.source)
    lab_path = f"{L}/consensus.csv" if a.source == "consensus" else f"{L}/candidate_labels.csv"
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
    X = T[feats].apply(pd.to_numeric, errors="coerce").values.astype(float); y = T.y.values; g = T.study_id.values
    # v5.2: leave-one-patient-out via css_model (imputation + scaling refitted inside every fold; L2 chosen by an inner
    # LOPO unless --no-nested; folds whose training part has one class are skipped and logged, never filled in)
    sc, fold_log = CM.lopo(X, y, g, a.l2, nested=not a.no_nested); T["score_lopo"] = sc
    base_col = "score_v4" if "score_v4" in T.columns else "score"
    # baselines on the same folds: the v4 rule score as it is, and one-family models (darkness only, vessel only)
    fam = {"darkness_only": [f for f in ("darkness_z", "edge_contrast") if f in feats],
           "vessel_only": [f for f in ("tube_ratio", "vessel_wm_mm", "local_cont_max_mm", "local_cont_min_mm",
                                       "terminal_ends", "net_radius_med_mm") if f in feats]}
    for k, fs in fam.items():
        T[f"score_{k}"] = CM.lopo(T[fs].apply(pd.to_numeric, errors="coerce").values.astype(float), y, g, a.l2)[0] if fs else np.nan
    ok = T.score_lopo.notna()
    fitted = {"score_lopo"} | {f"score_{k}" for k in fam}           # Brier only for fitted models, not rule scores
    def metrics(col):
        D = T[ok].assign(s=pd.to_numeric(T[col][ok], errors="coerce"))
        D = D[D.s.notna()]
        if not len(D): return dict(auc=None)
        r = lambda v: None if v is None or not np.isfinite(v) else round(float(v), 3)
        return dict(auc=r(CM.auc(D.y, D.s)), auc_ci95=CM.patient_bootstrap(D, lambda x: CM.auc(x.y, x.s), a.boot, a.seed),
                    pr_auc=r(CM.average_precision(D.y, D.s)), sens_at_1fp_per_scan=r(CM.sens_at_fp(D, "s", 1.0)),
                    sens_at_3fp_per_scan=r(CM.sens_at_fp(D, "s", 3.0)),
                    brier=r(CM.brier(D.y, D.s)) if col in fitted else None)
    evals = {"model": metrics("score_lopo"), base_col: metrics(base_col),
             **{k: metrics(f"score_{k}") for k, fs in fam.items() if fs}}
    auc_m, ci_m = evals["model"]["auc"], evals["model"].get("auc_ci95"); auc_b, ci_b = evals[base_col]["auc"], evals[base_col].get("auc_ci95")
    per = [dict(study_id=p, n=int((T.study_id == p).sum()), n_css=int(T.y[T.study_id == p].sum()),
                auc_model=CM.auc(T.y[T.study_id == p], T.score_lopo[T.study_id == p]),
                auc_score_v4=CM.auc(T.y[T.study_id == p], T[base_col][T.study_id == p])) for p in T.study_id.unique()]
    l2_final = a.l2
    if not a.no_nested:                                   # the final model's L2 by LOPO over all training studies
        best = (-1, a.l2)
        for cand in CM.L2_GRID:
            s_in = CM.lopo(X, y, g, cand)[0]; f_ = np.isfinite(s_in); v = CM.auc(y[f_], s_in[f_])
            if np.isfinite(v) and v > best[0]: best = (v, cand)
        l2_final = best[1]
    Z, med, mu, sd = prep(X); w, b = fit_lr(Z, y, l2_final)
    model = dict(model="logistic_l2_balanced", version="v5.2", fitted_at=pd.Timestamp.now().isoformat(timespec="seconds"),
                 features=feats, weights=[round(float(v), 6) for v in w], bias=round(b, 6),
                 impute_median=[float(v) for v in med], mean=[float(v) for v in mu], sd=[float(v) for v in sd],
                 l2=l2_final, l2_selection="nested LOPO over " + str(list(CM.L2_GRID)) if not a.no_nested else "fixed",
                 split=a.split, source=a.source, training_studies=sorted(T.study_id.unique().tolist()),
                 n_objects=int(len(T)), n_css=int(T.y.sum()), n_patients=int(n_pat),
                 lopo_auc=auc_m, lopo_auc_ci95=ci_m, baseline=base_col, baseline_auc=auc_b, baseline_auc_ci95=ci_b,
                 evaluation=evals, lopo_fold_status=fold_log,
                 per_patient=[{k: (round(v, 3) if isinstance(v, float) and np.isfinite(v) else (None if isinstance(v, float) else v))
                               for k, v in r.items()} for r in per],
                 dropped_features=gone,
                 # v5.1 provenance: detector builds of the training outputs, label file hash (rank_candidates checks)
                 detector_builds=sorted({str((json.load(open(f"{B}/review/{s}_manifest.json")).get("detect", {})
                                              .get("environment", {}) or {}).get("build", "unknown"))
                                         if os.path.exists(f"{B}/review/{s}_manifest.json") else "unknown"
                                         for s in T.study_id.unique()}),
                 label_file=lab_path, label_sha256=sha256(lab_path), build=build_id(), research_only=True,
                 note="research software; pooled LOPO AUC over few patients is noisy - read the intervals")
    model["model_hash"] = CM.model_hash({k: v for k, v in model.items() if k != "fitted_at"})
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(model, open(out, "w"), indent=2)
    top = sorted(zip(feats, w), key=lambda t: -abs(t[1]))[:8]
    print(f"LOPO AUC model {model['lopo_auc']} {ci_m} vs {base_col} {model['baseline_auc']} {ci_b} (patient bootstrap 95 %)")
    print("largest standardised weights: " + ", ".join(f"{f} {v:+.2f}" for f, v in top))
    for k, e in evals.items():
        print(f"  {k:14s} AUC {e.get('auc')} {e.get('auc_ci95')}  PR-AUC {e.get('pr_auc')}  sens@1FP/scan {e.get('sens_at_1fp_per_scan')}"
              + (f"  Brier {e['brier']}" if e.get("brier") is not None else ""))
    sk = [f["study_id"] for f in fold_log if f["status"] == "skipped"]
    if sk: print(f"  LOPO folds skipped (training part had one class): {', '.join(sk)}")
    card = write_card(model, os.path.join(os.path.dirname(out), "ranker_model_card.md"))
    print(f"model card -> {card} (model hash {model['model_hash'][:12]})")
    print(f"-> {out}   (apply: rank_candidates.py SUBJECT; never trained on test-split studies)")
