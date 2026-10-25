#!/usr/bin/env python3
"""Apply the trained ranking model to one study (v5.0, Layer 3).
usage: rank_candidates.py SUBJECT [--model review/ranker_model.json] [--allow-training-study]
Adds two columns at the END of review/ID_candidates.csv (replaced if present; nothing else changes - cand_id and
the existing score / score_v4 order stay, so both rankings are visible):
  score_model   model output 0-1 (NOT a calibrated probability of cSS)
  rank_model    1 = highest score_model
Features the model needs but this study lacks are filled with the training median (named in the output).
Refuses a study the model was trained on (its ranking would be optimistic) unless --allow-training-study.
CSS_QC_STRICT=1 and case QC low_confidence_review (review/ID_qc.json): score_model x 0.8 (a soft penalty that
keeps the order within the case; it matters only when candidates of several cases are pooled).
run_css.sh calls this only when the model exists and the study is not in the test split; a test-split study is
scored by running this script by hand with the frozen model. Research software, not for clinical use."""
import sys, os, json, argparse, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, FP, use_derived, sha256

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("subject"); ap.add_argument("--model", default="")
ap.add_argument("--allow-training-study", action="store_true")
a = ap.parse_args(); B = _base(); s = a.subject
mp = a.model or f"{B}/review/ranker_model.json"
if not os.path.exists(mp): sys.exit(f"no model ({mp}) - train one with train_ranker.py first")
M = json.load(open(mp))
if s in M.get("training_studies", []) and not a.allow_training_study:
    sys.exit(f"{s} was used to train {os.path.basename(mp)} - its model ranking would be optimistic "
             f"(--allow-training-study to rank anyway)")
# v5.1 stale-model checks: labels changed since the fit, or this study detected by a different detector build
strict = os.environ.get("CSS_QC_STRICT") == "1"
warn = []
lp = M.get("label_file")
if lp and os.path.exists(lp) and M.get("label_sha256") and sha256(lp) != M["label_sha256"]:
    warn.append(f"labels ({os.path.basename(lp)}) changed since the model was fitted - re-run train_ranker.py")
_mp = f"{B}/review/{s}_manifest.json"
_b = str((json.load(open(_mp)).get("detect", {}).get("environment", {}) or {}).get("build", "unknown")) if os.path.exists(_mp) else "unknown"
if M.get("detector_builds") and _b not in M["detector_builds"]:
    warn.append(f"{s} was detected with build {_b}, the model was trained on outputs of {', '.join(M['detector_builds'])}")
for w_ in warn: print(("ERROR: " if strict else "warning: ") + w_)
if warn and strict: sys.exit("CSS_QC_STRICT=1: stale model - not ranked")
cp = f"{B}/review/{s}_candidates.csv"
if not os.path.exists(cp): sys.exit(f"no {cp} - run detect_css.py {s} first")
d = pd.read_csv(cp).drop(columns=["score_model", "rank_model"], errors="ignore")
if not len(d):
    print(f"{s}: no candidates - nothing to rank"); d.assign(score_model=[], rank_model=[]).to_csv(cp, index=False); sys.exit(0)
X = pd.DataFrame(index=d.index)
sp = f"{B}/review/{s}_surface_features.csv"
sf = pd.read_csv(sp).drop(columns=["list", "obj_id"], errors="ignore").drop_duplicates(FP) \
    if use_derived(sp, [f"{B}/review/{s}_candidates.nii.gz"], "surface features") else None
full = d.merge(sf, on=FP, how="left", suffixes=("", "_surf")) if sf is not None else d
missing = [f for f in M["features"] if f not in full.columns]
for f in M["features"]: X[f] = pd.to_numeric(full[f], errors="coerce").values if f in full.columns else np.nan
Xv = X.values.astype(float); med = np.array(M["impute_median"])
Xv = np.where(np.isfinite(Xv), Xv, med)
t = ((Xv - np.array(M["mean"])) / np.array(M["sd"])) @ np.array(M["weights"]) + M["bias"]
sc = 1 / (1 + np.exp(-t))
note = ""
qp = f"{B}/review/{s}_qc.json"
if os.environ.get("CSS_QC_STRICT") == "1" and os.path.exists(qp) and json.load(open(qp)).get("status") == "low_confidence_review":
    sc = sc * 0.8; note = "  [QC low_confidence_review: score_model x0.8, CSS_QC_STRICT=1]"
d["score_model"] = np.round(sc, 4)
d["rank_model"] = pd.Series(sc).rank(ascending=False, method="first").astype(int).values
d.to_csv(cp, index=False)
order = d.sort_values("rank_model").cand_id.astype(int).tolist()
print(f"{s}: model ranking added (score_model, rank_model) - model fitted {M.get('fitted_at')} on "
      f"{M.get('n_patients')} patients, LOPO AUC {M.get('lopo_auc')} vs {M.get('baseline')} {M.get('baseline_auc')}{note}")
print(f"   order by model (cand_id): {order[:15]}{' ...' if len(order) > 15 else ''}")
if missing: print(f"   filled with the training median (not in this study's output): {', '.join(missing)}")
