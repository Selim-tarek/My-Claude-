#!/usr/bin/env python3
"""Score detector output against synthetic truth.  usage: eval_synthetic.py P005S
Reports sensitivity, worst rank, how many NON-lesion candidates (real veins/artifacts
of the host scan) rank above the worst lesion, and a ranking AUC."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
subj = sys.argv[1]
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
truth = nib.load(f"{base}/work/{subj}_truth.nii.gz").get_fdata().astype(int)
cand  = nib.load(f"{base}/review/{subj}_candidates.nii.gz").get_fdata().astype(int)
n = int(truth.max()); best, lesion_c = [], set()
for i in range(1, n + 1):
    ids = sorted(set(np.unique(cand[ndi.binary_dilation(truth == i, iterations=1)])) - {0})
    if ids: best.append(ids[0]); lesion_c |= set(ids)
    print(f"  lesion {i}: {'FOUND  (candidate #' + str(ids[0]) + ')' if ids else 'MISSED'}")
all_c = set(np.unique(cand)) - {0}
non = sorted(all_c - lesion_c)
worst = max(best) if best else 0
above = sum(1 for c in non if c < worst)
pairs = [(b < c) for b in best for c in non]
auc = float(np.mean(pairs)) if pairs else (1.0 if best else 0.0)
print(f"\n  sensitivity: {len(best)}/{n} = {100*len(best)/max(n,1):.0f}%   worst rank #{worst}")
print(f"  non-lesion candidates ranked above the worst lesion: {above}   ranking AUC: {auc:.2f}")
print(f"RESULT found={len(best)} total={n} worst={worst} above={above} auc={auc:.3f} ncand={len(all_c)}")
