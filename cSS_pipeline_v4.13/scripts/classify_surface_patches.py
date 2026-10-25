#!/usr/bin/env python3
"""cSS-versus-vein ranking from surface-aligned patches (v5.2, master plan Stage 3). Research software: the output
ranks candidates for expert review; it is not a diagnosis and not a calibrated probability.
usage:
  classify_surface_patches.py features SUBJECT [SUBJECT ...]
        review/ID_patches/*.npz (surface_aligned_patch.py) -> review/ID_patch_features.csv
  classify_surface_patches.py train [--split dev] [--source consensus|reader:NAME] [--negatives vein|all]
        [--min-studies 5] [--l2 1.0] [--no-nested] [--out review/surface_patch_model.json]
  classify_surface_patches.py predict SUBJECT [SUBJECT ...] [--model review/surface_patch_model.json]
        [--allow-training-study]          -> review/ID_patch_predictions.csv
  classify_surface_patches.py lopo-manifest MANIFEST.csv [--l2 1.0]
        explicit patch list (columns patch_path, study_id, label = cSS | vein code): LOPO by study, no model saved
Library: read_manifest / build_matrix / grouped_lopo (same features, same css_model LOPO as 'train').
Features per channel, in a FIXED order (swi, phase, flair, labels): validity fraction, robust-z distribution
(mean, sd, p5, p50, p95, min, max), mean |gradient| along u / v / n, and 8-bin profiles (mean, sd) along each axis;
for the label channel the tissue share in 8 bins along the normal. Every channel adds <ch>_missing and
<ch>_valid_fraction: a missing phase / FLAIR is information and is never a zero image. Optional geometry and vessel
columns are read from review/ID_candidates|excluded.csv (surface / vessel topology, v5.1.x).
Training: labelled objects of the DEV split only (labels/study_meta.csv; --split test is refused), target cSS (1)
vs vein (VS / VC / VT / VX / V; --negatives all = every non-cSS label), 'uncertain' left out. Class-balanced L2
logistic regression, leave-one-patient-out with nested L2 selection (css_model.py); folds whose training part has
one class are skipped and reported. Ablations on the same folds: A swi, B swi+phase, C swi+labels+surface geometry,
D swi+vessel topology, E swi+flair, F swi+phase+labels+geometry+vessel, G everything. Reported per ablation: AUC,
PR-AUC, Brier, sensitivity at 1 FP per scan, with patient-bootstrap intervals for AUC. Needs >= --min-studies studies
and both classes, else it refuses. The model JSON records the feature schema, channel order, preprocessing,
training studies, fold log, out-of-fold predictions, label-file hash, build and a model hash."""
import sys, os, json, argparse, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, label_class, sha256, build_id
import css_model as CM

FEATURE_SCHEMA = "css-patch-features/1.0"
CHANNELS = ["swi", "phase", "flair", "labels"]
GEOMETRY = ["pial_dist_mm", "csf_depth_mm", "sulcal_depth_mm", "surface_contact_frac_1mm", "elongation_3d",
            "compactness_3d", "surf_width_mm", "surf_inradius_mm", "surf_aspect", "geo_length_mm"]
VESSEL = ["tube_ratio", "local_cont_max_mm", "local_cont_min_mm", "terminal_ends", "local_branch_per10mm",
          "vessel_wm_mm", "vessel_run_mm", "tortuosity", "cand_radius_med_mm", "net_radius_med_mm"]
ABLATIONS = {"A_swi": ["swi"], "B_swi_phase": ["swi", "phase"], "C_swi_surface": ["swi", "labels", "geometry"],
             "D_swi_vessel": ["swi", "vessel"], "E_swi_flair": ["swi", "flair"],
             "F_swi_phase_surface_vessel": ["swi", "phase", "labels", "geometry", "vessel"],
             "G_all": ["swi", "phase", "flair", "labels", "geometry", "vessel"]}


def robust_z(x, valid):
    valid = valid & np.isfinite(x); z = np.full(x.shape, np.nan, np.float32)
    if not valid.any(): return z
    med = np.median(x[valid]); mad = np.median(np.abs(x[valid] - med))
    z[valid] = (x[valid] - med) / max(1.4826 * mad, float(np.std(x[valid])), 1e-3)
    return z


def axis_profile(z, valid, axis, bins=8):
    e = np.linspace(0, z.shape[axis], bins + 1).astype(int); m, s = [], []
    for lo, hi in zip(e[:-1], e[1:]):
        sl = [slice(None)] * 3; sl[axis] = slice(lo, max(lo + 1, hi))
        v = z[tuple(sl)][valid[tuple(sl)]]
        m.append(float(v.mean()) if len(v) else np.nan); s.append(float(v.std()) if len(v) else np.nan)
    return m + s


def channel_features(name, img, valid):
    """-> ordered (names, values); img None = channel missing (values NaN, <ch>_missing = 1)"""
    names = [f"{name}_missing", f"{name}_valid_fraction"]
    if name == "labels":
        names += [f"labels_tissue_n{b}" for b in range(8)]
        if img is None: return names, [1.0, 0.0] + [np.nan] * 8
        t = ((img > 0) & (img != 24)).astype(np.float32)
        return names, [0.0, float(valid.mean())] + axis_profile(np.where(valid, t, np.nan), valid, 2)[:8]
    stats = ["mean", "sd", "p5", "p50", "p95", "min", "max", "grad_u", "grad_v", "grad_n"]
    prof = [f"prof_{ax}_{k}{b}" for ax in "uvn" for k in ("mean", "sd") for b in range(8)]
    names += [f"{name}_{s}" for s in stats + prof]
    if img is None or not valid.any(): return names, [1.0 if img is None else 0.0, float(valid.mean()) if img is not None else 0.0] + [np.nan] * (len(names) - 2)
    z = robust_z(img.astype(np.float32), valid); v = z[valid]
    g = [float(np.nanmean(np.abs(gg[valid]))) for gg in np.gradient(np.where(valid, z, 0.0))]
    vals = [0.0, float(valid.mean()), float(v.mean()), float(v.std()), *map(float, np.percentile(v, [5, 50, 95])),
            float(v.min()), float(v.max()), *g]
    for ax in range(3): vals += axis_profile(z, valid, ax)
    return names, vals


def patch_features(npz_path):
    z = np.load(npz_path, allow_pickle=False); names, vals = [], []
    for ch in CHANNELS:
        img = z[f"image_{ch}"] if f"image_{ch}" in z.files else None
        vm = z[f"valid_{ch}"] if f"valid_{ch}" in z.files else np.zeros(z["valid"].shape, bool)
        n_, v_ = channel_features(ch, img, np.asarray(vm).astype(bool))
        names += n_; vals += v_
    return names, vals, json.loads(str(z["metadata_json"])) if "metadata_json" in z.files else {}


POSITIVE = {"cSS", "css", "C"}
NEGATIVE = {"vein", "surface_vein", "other_vein", "V", "VS", "VC", "VT", "VX"}


def read_manifest(path):
    """explicit patch list -> rows (dicts). Needs columns patch_path, study_id, label; a relative patch_path is read
    against the manifest's folder, then $CSS_BASE. A patch listed twice (whatever its labels) or a missing file -> ValueError."""
    d = pd.read_csv(path, dtype=str)
    need = {"patch_path", "study_id", "label"} - set(d.columns)
    if need: raise ValueError(f"{path}: missing columns {sorted(need)}")
    here = os.path.dirname(os.path.abspath(path)); rows, seen = [], set()
    for _, r in d.iterrows():
        pp = r.patch_path
        if not os.path.isabs(pp):
            pp = next((c for c in (os.path.join(here, pp), os.path.join(_base(), pp)) if os.path.exists(c)), os.path.join(here, pp))
        key = os.path.realpath(pp)
        if key in seen: raise ValueError(f"{path}: patch listed more than once: {r.patch_path}")
        if not os.path.exists(pp): raise ValueError(f"{path}: patch file not found: {r.patch_path}")
        seen.add(key); rows.append(dict(r, patch_path=pp))
    return rows


def build_matrix(rows):
    """rows from read_manifest -> (X, y, groups, feature_names, channels). Features of channels present in at least one
    patch, named '<channel>:<feature>' (e.g. 'phase:missing' = 1 where that patch has no phase). Labels: cSS -> 1,
    vein codes -> 0; anything else -> ValueError (no silent relabelling)."""
    bad = sorted({r["label"] for r in rows if r["label"] not in POSITIVE | NEGATIVE})
    if bad: raise ValueError(f"labels must be cSS or a vein code, got {bad}")
    feats, present = [], set()
    for r in rows:
        names, vals, _ = patch_features(r["patch_path"]); feats.append(vals)
        z = np.load(r["patch_path"], allow_pickle=False); present |= {ch for ch in CHANNELS if f"image_{ch}" in z.files}
    keep = [k for k, n in enumerate(names) if n.split("_", 1)[0] in present]
    X = np.asarray(feats, float)[:, keep]
    fn = [names[k].replace("_", ":", 1) for k in keep]
    y = np.array([1 if r["label"] in POSITIVE else 0 for r in rows]); g = np.array([r["study_id"] for r in rows])
    return X, y, g, fn, sorted(present)


def grouped_lopo(X, y, groups, l2=1.0, nested=False):
    """leave-one-study-out scores + fold log (css_model.lopo: folds by study, one-class training folds skipped)"""
    return CM.lopo(X, y, groups, l2, nested=nested)


auc = CM.auc


def cmd_features(subjects, B):
    for s in subjects:
        mp = f"{B}/review/{s}_patches.csv"
        if not os.path.exists(mp): print(f"{s}: no {mp} - run surface_aligned_patch.py first"); continue
        man = pd.read_csv(mp); rows = []
        for _, r in man[man.status == "ok"].iterrows():
            names, vals, meta = patch_features(f"{B}/{r.patch_path}")
            rows.append(dict(study_id=s, cand_uid=r.cand_uid, list=r.list, obj_id=r.obj_id, **dict(zip(names, vals)),
                             feature_schema=FEATURE_SCHEMA))
        pd.DataFrame(rows).to_csv(f"{B}/review/{s}_patch_features.csv", index=False)
        print(f"{s}: {len(rows)} patch feature rows ({len(names) if rows else 0} features) -> review/{s}_patch_features.csv")


def study_table(B, s):
    p = f"{B}/review/{s}_patch_features.csv"
    if not os.path.exists(p): return None
    d = pd.read_csv(p)
    geo = []
    for lst in ("candidates", "excluded"):
        q = f"{B}/review/{s}_{lst}.csv"
        if os.path.exists(q):
            g = pd.read_csv(q); g["list"] = lst; g["obj_id"] = g["cand_id" if lst == "candidates" else "excl_id"]; geo.append(g)
    if geo:
        g = pd.concat(geo, ignore_index=True)
        sp = f"{B}/review/{s}_surface_features.csv"
        if os.path.exists(sp):
            sf = pd.read_csv(sp); sf["list"] = sf["list"].map({"candidate": "candidates", "excluded": "excluded"})
            g = g.merge(sf[["list", "obj_id"] + [c for c in GEOMETRY if c in sf.columns and c not in g.columns]],
                        on=["list", "obj_id"], how="left")
        keep = ["list", "obj_id"] + [c for c in GEOMETRY + VESSEL if c in g.columns]
        d = d.merge(g[keep], on=["list", "obj_id"], how="left")
    return d


def columns_for(groups, cols):
    out = []
    for gname in groups:
        if gname in CHANNELS: out += [c for c in cols if c.startswith(gname + "_")]
        elif gname == "geometry": out += [c for c in GEOMETRY if c in cols]
        elif gname == "vessel": out += [c for c in VESSEL if c in cols]
    return out


def cmd_train(a, B):
    if a.split == "test": sys.exit("training on the test split is refused (labels/study_meta.csv): test studies are for reporting")
    L = f"{B}/labels"; meta = pd.read_csv(f"{L}/study_meta.csv", dtype=str)
    studies = meta[meta.split == "dev"].study_id if a.split == "dev" else meta[meta.split != "test"].study_id
    if a.source == "consensus":
        lab = pd.read_csv(f"{L}/consensus.csv"); lab = lab[lab.method != "pending"].assign(cls=lambda x: x.consensus_class)
        lab_path = f"{L}/consensus.csv"
    else:
        rd = a.source.split(":", 1)[1]; lab = pd.read_csv(f"{L}/candidate_labels.csv"); lab = lab[lab.reader_id == rd]
        lab = lab[lab.session == lab.groupby("study_id").session.transform("max")].assign(cls=lambda x: x.label.map(label_class))
        lab_path = f"{L}/candidate_labels.csv"
    neg = {"surface_vein", "other_vein"} if a.negatives == "vein" else None
    rows = []
    for s in studies:
        d = study_table(B, s)
        if d is None: continue
        m = d.merge(lab[["study_id", "cand_uid", "cls"]], on=["study_id", "cand_uid"], how="inner")
        m = m[(m.cls == "cSS") | (m.cls.isin(neg) if neg else (m.cls != "uncertain"))]
        if len(m): rows.append(m.assign(y=(m.cls == "cSS").astype(int)))
    if not rows: sys.exit("no labelled patches in the requested split - extract patches and label studies first")
    T = pd.concat(rows, ignore_index=True); n_st = T.study_id.nunique()
    print(f"train: {len(T)} labelled patches ({int(T.y.sum())} cSS, {int((1 - T.y).sum())} "
          f"{'vein' if neg else 'non-cSS'}) from {n_st} studies, split {a.split}")
    if n_st < a.min_studies: sys.exit(f"need >= {a.min_studies} labelled studies, have {n_st} - not fitted")
    if T.y.nunique() < 2: sys.exit("need both cSS and negative labels - not fitted")
    cols = [c for c in T.columns if c not in ("study_id", "cand_uid", "list", "obj_id", "cls", "y", "feature_schema")]
    res, oof = {}, {}
    for name, groups in ABLATIONS.items():
        fc = [c for c in columns_for(groups, cols) if T[c].notna().any()]
        if not fc: res[name] = dict(n_features=0, note="no features available"); continue
        sc, log = CM.lopo(T[fc].values.astype(float), T.y.values, T.study_id.values, a.l2, nested=not a.no_nested)
        ok = np.isfinite(sc); D = T[ok].assign(s=sc[ok])
        res[name] = dict(n_features=len(fc), auc=CM.auc(D.y, D.s), auc_ci95=CM.patient_bootstrap(D, lambda x: CM.auc(x.y, x.s), 300),
                         pr_auc=CM.average_precision(D.y, D.s), brier=CM.brier(D.y, D.s),
                         sens_at_1fp_per_scan=CM.sens_at_fp(D, "s", 1.0),
                         folds_skipped=[f["study_id"] for f in log if f["status"] == "skipped"])
        oof[name] = (sc, log)
    fc = [c for c in columns_for(ABLATIONS["G_all"], cols) if T[c].notna().any()]
    Z, med, mu, sd = CM.prep(T[fc].values.astype(float)); w, b = CM.fit_lr(Z, T.y.values, a.l2)
    model = dict(model_type="logistic_l2_balanced", model_version="v5.2", feature_schema=FEATURE_SCHEMA, channel_order=CHANNELS,
                 feature_names=fc, weights=[float(x) for x in w], bias=b, imputation_medians=[float(x) for x in med],
                 training_mean=[float(x) for x in mu], training_scale=[float(x) for x in sd], l2=a.l2,
                 training_studies=sorted(T.study_id.unique().tolist()), n_patches=int(len(T)), n_css=int(T.y.sum()),
                 n_negative=int((1 - T.y).sum()), negatives=a.negatives, split=a.split, source=a.source,
                 label_file=lab_path, label_sha256=sha256(lab_path), ablations=res,
                 lopo_fold_status=oof["G_all"][1] if "G_all" in oof else [],
                 lopo_predictions=[dict(study_id=s_, cand_uid=u_, y=int(y_), score=None if not np.isfinite(p_) else round(float(p_), 4))
                                   for s_, u_, y_, p_ in zip(T.study_id, T.cand_uid, T.y, oof["G_all"][0])] if "G_all" in oof else [],
                 build=build_id(), research_only=True,
                 note="ranks candidates for expert review; not calibrated, not clinically validated")
    model["model_hash"] = CM.model_hash(model)
    out = a.out or f"{B}/review/surface_patch_model.json"; json.dump(model, open(out, "w"), indent=2, default=float)
    for name, r in res.items():
        if "auc" in r:
            print(f"  {name:28s} n_feat {r['n_features']:4d}  AUC {r['auc']:.3f} {r['auc_ci95']}  PR-AUC {r['pr_auc']:.3f}  "
                  f"Brier {r['brier']:.3f}  sens@1FP {r['sens_at_1fp_per_scan']:.2f}" + (f"  skipped {r['folds_skipped']}" if r["folds_skipped"] else ""))
    print(f"-> {out} (model hash {model['model_hash'][:12]}; fitted on {len(model['training_studies'])} dev studies)")


def cmd_predict(a, B):
    mp = a.model or f"{B}/review/surface_patch_model.json"
    if not os.path.exists(mp): sys.exit(f"no model {mp}")
    M = json.load(open(mp))
    if M.get("feature_schema") != FEATURE_SCHEMA or M.get("channel_order") != CHANNELS:
        sys.exit(f"model feature schema / channel order ({M.get('feature_schema')}, {M.get('channel_order')}) do not match this "
                 f"script ({FEATURE_SCHEMA}, {CHANNELS}) - retrain")
    if M.get("model_hash") != CM.model_hash(M): sys.exit("model hash mismatch - the model file was changed after fitting")
    for s in a.subjects:
        if s in M["training_studies"] and not a.allow_training_study:
            print(f"{s}: used to train the model - not predicted (--allow-training-study)"); continue
        d = study_table(B, s)
        if d is None: print(f"{s}: no patch features - run 'features' first"); continue
        miss = [c for c in M["feature_names"] if c not in d.columns]
        if [c for c in miss if not any(c.startswith(ch + "_") for ch in CHANNELS)]:
            sys.exit(f"{s}: features the model needs are absent: {miss[:5]} - schema mismatch")
        X = np.column_stack([pd.to_numeric(d[c], errors="coerce").values if c in d.columns else np.full(len(d), np.nan)
                             for c in M["feature_names"]])
        sc = CM.predict(X, np.array(M["imputation_medians"]), np.array(M["training_mean"]), np.array(M["training_scale"]),
                        np.array(M["weights"]), M["bias"])
        mm = d[[c for c in d.columns if c.endswith("_missing")]]
        out = pd.DataFrame(dict(study_id=s, cand_uid=d.cand_uid, list=d.list, obj_id=d.obj_id, css_score=np.round(sc, 4),
                                missing_modalities=[";".join(c[:-8] for c in mm.columns if r_[c] == 1) for _, r_ in mm.iterrows()],
                                model_hash=M["model_hash"], model_version=M["model_version"], build=build_id(), research_only=1))
        out.to_csv(f"{B}/review/{s}_patch_predictions.csv", index=False)
        print(f"{s}: {len(out)} patch scores -> review/{s}_patch_predictions.csv (model {M['model_hash'][:12]}; research only)")


def cmd_lopo_manifest(a):
    X, y, g, fn, ch = build_matrix(read_manifest(a.manifest))
    sc, log = grouped_lopo(X, y, g, a.l2); ok = np.isfinite(sc)
    sk = [f["study_id"] for f in log if f["status"] == "skipped"]
    print(f"{len(y)} patches ({int(y.sum())} cSS) from {len(set(g))} studies, channels {ch}, {len(fn)} features: LOPO AUC "
          f"{CM.auc(y[ok], sc[ok]):.3f}, PR-AUC {CM.average_precision(y[ok], sc[ok]):.3f}" + (f"; folds skipped {sk}" if sk else "")
          + "  [research only, no model saved]")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("features"); f.add_argument("subjects", nargs="+")
    t = sub.add_parser("train"); t.add_argument("--split", default="dev", choices=["dev", "all", "test"])
    t.add_argument("--source", default="consensus"); t.add_argument("--negatives", default="vein", choices=["vein", "all"])
    t.add_argument("--min-studies", type=int, default=5); t.add_argument("--l2", type=float, default=1.0)
    t.add_argument("--no-nested", action="store_true"); t.add_argument("--out", default="")
    p = sub.add_parser("predict"); p.add_argument("subjects", nargs="+"); p.add_argument("--model", default="")
    p.add_argument("--allow-training-study", action="store_true")
    q = sub.add_parser("lopo-manifest"); q.add_argument("manifest"); q.add_argument("--l2", type=float, default=1.0)
    a = ap.parse_args(); B = _base()
    {"features": lambda: cmd_features(a.subjects, B), "train": lambda: cmd_train(a, B), "predict": lambda: cmd_predict(a, B),
     "lopo-manifest": lambda: cmd_lopo_manifest(a)}[a.cmd]()
