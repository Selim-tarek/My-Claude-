#!/usr/bin/env python3
"""Safety check: are KNOWN lesions still detected after a detector change?
usage: check_known.py SUBJECT --old OLD_candidates.nii.gz --ids 3,4,5
For each listed lesion of an older candidate map (e.g. review/v3/P006_candidates.nii.gz) it reports
the current candidate(s) covering it, or the exclusion reason (review/ID_excluded.csv), or
"NOT DETECTED". Grids must match (same imported SWI)."""
import sys, os, numpy as np, pandas as pd, nibabel as nib
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base

subj = sys.argv[1]; base = _base()
arg = lambda f: sys.argv[sys.argv.index(f) + 1] if f in sys.argv else ""
if not arg("--old") or not arg("--ids"): sys.exit(__doc__)
old = nib.load(os.path.expanduser(arg("--old"))).get_fdata().astype(int)
cand = nib.load(f"{base}/review/{subj}_candidates.nii.gz").get_fdata().astype(int)
xp = f"{base}/review/{subj}_excluded.nii.gz"
excl = nib.load(xp).get_fdata().astype(int) if os.path.exists(xp) else np.zeros_like(cand)
xr = pd.read_csv(f"{base}/review/{subj}_excluded.csv").set_index("excl_id")["reason"] \
    if os.path.exists(f"{base}/review/{subj}_excluded.csv") else pd.Series(dtype=str)
if old.shape != cand.shape: sys.exit("grids differ - re-import the same SWI")
kept = 0
for i in [int(x) for x in arg("--ids").split(",") if x.strip()]:
    m = ndi.binary_dilation(old == i, iterations=2)
    c = sorted(set(np.unique(cand[m])) - {0}); x = sorted(set(np.unique(excl[m])) - {0})
    if c: kept += 1
    msg = (f"kept as candidate {', '.join('#' + str(v) for v in c)}" if c else "") + \
          ("; " if c and x else "") + ("; ".join(f"EXCLUDED ({xr.get(v, '?')})" for v in x) if x else "")
    print(f"known lesion #{i}: {msg or 'NOT DETECTED'}")
print(f"{subj}: {kept} of {len(arg('--ids').split(','))} known lesions still shown for review")
