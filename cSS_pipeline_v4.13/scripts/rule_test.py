#!/usr/bin/env python3
"""Which single-feature exclusion rule would remove false positives WITHOUT losing known cSS?
usage: rule_test.py --pos P006 --old OLD_candidates.nii.gz --ids 3,4,5,13,18,19 --neg P011[,P010]
                    [--synthetic results/stress_TAG_features.csv] [--top 12]
  --pos/--old/--ids  known cSS: lesions of an older candidate map (as in check_known.py), matched to
                     the CURRENT candidates of that subject
  --neg              cSS-negative subjects: every current candidate there is a false positive
  --synthetic        optional stress-test feature table: rules are also checked on synthetic lesions
                     (is_lesion == 1); synthetic features are partly artifacts of the generator
For each feature and direction the threshold is put just beyond the most extreme known lesion (plus a
safety margin), so every known lesion is kept; the table shows how many false positives that removes.
With only a few known lesions this is a SHORTLIST for discussion, not a validated rule."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, numpy as np, pandas as pd, nibabel as nib
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base

base = _base()
arg = lambda f, d="": sys.argv[sys.argv.index(f) + 1] if f in sys.argv else d
if not (arg("--pos") and arg("--old") and arg("--ids") and arg("--neg")): sys.exit(__doc__)
pos_s, top = arg("--pos"), int(arg("--top", "12"))
FEATS = ["darkness_z", "volume_mm3", "elongation", "surface_contact", "cortex_frac", "cortex_dist_mm",
         "surface_gradient", "pial_dist_mm", "bank_frac", "tube_ratio", "surface_alignment", "tram_frac",
         "vein_tree_mm", "mirror_dark_frac", "edge_contrast", "parenchyma_frac", "extent_mm",
         "vessel_ext_mm", "vessel_wm_mm", "vessel_run_mm", "axis_normal", "css_evidence", "vein_evidence",
         "score_v4"]

# known lesions -> current candidates of the positive subject
old = nib.load(os.path.expanduser(arg("--old"))).get_fdata().astype(int)
cand = nib.load(f"{base}/review/{pos_s}_candidates.nii.gz").get_fdata().astype(int)
if old.shape != cand.shape: sys.exit("grids differ - re-import the same SWI")
dp = pd.read_csv(f"{base}/review/{pos_s}_candidates.csv")
known = {}
for i in [int(x) for x in arg("--ids").split(",") if x.strip()]:
    v = cand[ndi.binary_dilation(old == i, iterations=2)]; v = v[v > 0]
    if len(v): known[i] = int(np.bincount(v).argmax())          # the candidate covering it most
    else: print(f"known lesion #{i}: not a current candidate (excluded or not detected) - not used")
if not known: sys.exit("no known lesion is a current candidate")
P = dp[dp.cand_id.isin(set(known.values()))]
O = dp[~dp.cand_id.isin(set(known.values()))]          # other candidates of the positive subject
N = pd.concat([pd.read_csv(f"{base}/review/{s}_candidates.csv").assign(subject=s)
               for s in arg("--neg").split(",")], ignore_index=True)
S = None
if arg("--synthetic"):
    S = pd.read_csv(os.path.expanduser(arg("--synthetic"))); S = S[S.is_lesion == 1]
F = [f for f in FEATS if f in dp.columns and f in N.columns]

print(f"known cSS ({pos_s}): {len(P)} candidates {sorted(known.items())}  (old id: current id)")
print(f"false positives: {len(N)} candidates of {arg('--neg')}" +
      (f";  synthetic lesions: {len(S)}" if S is not None else "") +
      f";  other {pos_s} candidates (unlabelled): {len(O)}\n")
print("Feature values of the known lesions vs the false positives (median [min-max]):")
for f in F:
    p_, n_ = P[f].dropna(), N[f].dropna()
    if not len(p_) or not len(n_): continue
    print(f"  {f:18s} known {p_.median():7.2f} [{p_.min():7.2f} {p_.max():7.2f}]   "
          f"FP {n_.median():7.2f} [{n_.min():7.2f} {n_.max():7.2f}]")

rows = []
for f in F:
    p_ = P[f].dropna()
    if len(p_) < len(P) or not N[f].notna().any(): continue      # a rule must be defined for every known lesion
    allv = pd.concat([p_, N[f].dropna()]); margin = 0.1 * (allv.max() - allv.min())
    for op, t in ((">", p_.max() + margin), ("<", p_.min() - margin)):
        hit = (lambda x: x > t) if op == ">" else (lambda x: x < t)
        nrem = int(hit(N[f].dropna()).sum())
        if nrem == 0: continue
        rows.append(dict(rule=f"exclude if {f} {op} {t:.2f}", fp_removed=nrem,
                         fp_pct=round(100 * nrem / len(N)),
                         other_pos_removed=int(hit(O[f].dropna()).sum()),
                         synthetic_lost=(f"{int(hit(S[f].dropna()).sum())}/{len(S)}" if S is not None and f in S else "-"),
                         margin_to_known=round(margin, 2)))
print("\nRules that keep ALL known lesions (threshold = most extreme known lesion + 10 % of the range):")
if rows:
    R = pd.DataFrame(rows).sort_values("fp_removed", ascending=False).head(top)
    print(R.to_string(index=False))
else:
    print("  none - every feature of the false positives overlaps the known lesions")
print("\nRead with care: few known lesions -> a rule can look safe by chance. Prefer rules that also lose "
      "few synthetic lesions and that make sense radiologically; confirm on the next patients.")
