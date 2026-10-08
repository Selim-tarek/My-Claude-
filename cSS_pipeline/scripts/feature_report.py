#!/usr/bin/env python3
"""Compare candidate features between calls (cSS vs Vein/Normal/Artifact) across all reviewed subjects.
usage: feature_report.py [--source expert|reader]
Default: for each subject the EXPERT reference (review/ID_expert_*.csv from expert_import.py) is used
when it exists, otherwise the reader's own calls (review/ID_calls.csv from review_css.py)."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, load_calls, match_calls
base = _base()
src_pref = sys.argv[sys.argv.index("--source") + 1] if "--source" in sys.argv else "expert"
subjects = sorted({os.path.basename(f).split("_calls.csv")[0] for f in glob.glob(f"{base}/review/*_calls.csv")} |
                  {os.path.basename(f).split("_expert_")[0] for f in glob.glob(f"{base}/review/*_expert_*.csv")
                   if not f.endswith("_expert_key.csv")})
rows = []; used = {}
for s in subjects:
    cp = f"{base}/review/{s}_candidates.csv"
    if not os.path.exists(cp): continue
    d = pd.read_csv(cp); calls = {}
    ex = sorted(f for f in glob.glob(f"{base}/review/{s}_expert_*.csv") if not f.endswith("_expert_key.csv"))
    if ex and src_pref == "expert":
        calls, dropped = match_calls(d, pd.read_csv(ex[0])); calls = calls or {}
        used[s] = "expert " + os.path.basename(ex[0]).split("_expert_")[1][:-4]
        if dropped: print(f"WARNING: {s}: {dropped} expert calls match no current candidate")
    else:
        calls, warn = load_calls(base, s, d)          # matched by fingerprint, not by rank
        if warn: print("WARNING:", warn)
        used[s] = "reader"
    if not calls: continue
    m = d.merge(pd.DataFrame(list(calls.items()), columns=["cand_id", "call"]), on="cand_id")
    m["subject"] = s; rows.append(m)
print("labels: " + ", ".join(f"{s}={u}" for s, u in used.items()))
if not rows: sys.exit("no reviewed cases yet (run review_css.py first)")
D = pd.concat(rows, ignore_index=True)
D = D[D.call.isin(["cSS", "Vein", "Normal", "Artifact"])]
print(f"{D.subject.nunique()} subject(s), {len(D)} decided candidates: {D.call.value_counts().to_dict()}\n")
def auc(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float); x, y = x[~np.isnan(x)], y[~np.isnan(y)]
    if not len(x) or not len(y): return float("nan")
    return float((x[:, None] > y[None, :]).mean() + 0.5 * (x[:, None] == y[None, :]).mean())
feats = ["darkness_z", "volume_mm3", "n_slices", "elongation", "surface_contact", "cortex_frac",
         "cortex_dist_mm", "branch_per10mm", "surface_gradient",
         "pial_dist_mm", "bank_frac", "tube_ratio", "surface_alignment", "tram_frac", "vein_tree_mm",
         "mirror_dark_frac", "flair_csf_z", "parenchyma_frac", "css_evidence", "vein_evidence",
         "score_v3", "score_v4", "score"]
print(f"{'feature':18s} {'cSS':>8s} {'Vein':>8s} {'Normal':>8s} {'Artifact':>9s}   AUC cSS vs rest")
for f in feats:
    if f not in D: continue
    med = D.groupby("call")[f].median()
    a = auc(D[D.call == "cSS"][f], D[D.call != "cSS"][f])
    tag = "  <- useful" if a >= 0.75 or a <= 0.25 else ""
    print(f"{f:18s} " + " ".join(f"{med.get(k, np.nan):8.2f}" for k in ["cSS", "Vein", "Normal"]) +
          f" {med.get('Artifact', np.nan):9.2f}   {a:5.2f}{tag}")
print("\nAUC near 0.5 = no help; >0.75 = higher in cSS; <0.25 = lower in cSS.")
