#!/usr/bin/env python3
"""Blinded picture sheet of all candidates for an expert reader -> review/ID_expert_sheet.pdf
usage: expert_sheet.py SUBJECT [--top N] [--seed 1] [--include-excluded]
  --include-excluded (v4.15): the automatically EXCLUDED objects (review/ID_excluded.csv) are mixed in,
      unmarked, so the expert's calls show whether an exclusion rule removed real cSS and give
      labelled hard negatives (surface veins, artifacts, ...). Recommended for validation cases.

Each candidate gets a SHEET NUMBER in random order (the detector's rank is hidden, so the
order cannot bias the expert) and one row of pictures: SWI zoom (40 mm), the same with the
candidate outlined in red, and an 8 mm minIP slab (veins = continuous tubes). No features or
scores are printed. Page 1 explains the task; the last page is an answer grid.
The key (sheet number -> cand_id + fingerprint) is saved separately in
review/ID_expert_key.csv - do not give it to the expert.
Import the answers with expert_import.py."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, numpy as np, nibabel as nib, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, FP, LABEL_HELP

subj = sys.argv[1]
top = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else None
seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 1
base = _base()
swi = nib.as_closest_canonical(nib.load(f"{base}/data/{subj}_swi.nii"))
cnd = nib.as_closest_canonical(nib.load(f"{base}/review/{subj}_candidates.nii.gz"))
I = swi.get_fdata().astype(np.float32); C = cnd.get_fdata().astype(int)
vox = swi.header.get_zooms()[:3]
br = I > 0
lo, hi = np.percentile(I[br], [1, 99]) if br.any() else (I.min(), I.max())
df = pd.read_csv(f"{base}/review/{subj}_candidates.csv")
objs = [("candidate", int(c)) for c in df.cand_id.astype(int).tolist()[:top]]
X = None
if "--include-excluded" in sys.argv and os.path.exists(f"{base}/review/{subj}_excluded.csv"):
    xdf = pd.read_csv(f"{base}/review/{subj}_excluded.csv")
    if len(xdf):
        X = nib.as_closest_canonical(nib.load(f"{base}/review/{subj}_excluded.nii.gz")).get_fdata().astype(int)
        objs += [("excluded", int(e)) for e in xdf.excl_id]
order = np.random.default_rng(seed).permutation(len(objs))
sheet = {n + 1: objs[k] for n, k in enumerate(order)}          # sheet number -> (list, id)
rows = []
for n in sorted(sheet):
    lst, i = sheet[n]
    r = (df.set_index("cand_id") if lst == "candidate" else xdf.set_index("excl_id")).loc[i]
    rows.append(dict(sheet_no=n, cand_id=i if lst == "candidate" else None, **{f: r[f] for f in FP},
                     list=lst, obj_id=i))
key = pd.DataFrame(rows, columns=["sheet_no", "cand_id"] + FP + ["list", "obj_id"])
key.to_csv(f"{base}/review/{subj}_expert_key.csv", index=False)
mask_of = lambda n: (C == sheet[n][1]) if sheet[n][0] == "candidate" else (X == sheet[n][1])

HW = int(round(20 / vox[0]))                                   # 40 mm window
SLAB = max(1, int(round(4 / vox[2])))
PER_PAGE = 4
out = f"{base}/review/{subj}_expert_sheet.pdf"
with PdfPages(out) as pdf:
    f = plt.figure(figsize=(8.27, 11.69)); f.text(0.08, 0.92, f"cSS expert reading - {subj}", fontsize=16, weight="bold")
    f.text(0.06, 0.86, (
        "For each numbered candidate write one code (optionally ':1'-':5' = confidence,\n"
        "5 = certain; e.g. 12VS:4):\n\n"
        + "".join(f"   {c:<22s}{t}\n" for c, t in LABEL_HELP) + "\n"
        "Single letters C / V / N / H / U are enough if you prefer.\n\n"
        "Definition used (Charidimou et al., Neurology 2017): well-defined, homogeneous\n"
        "hypointense curvilinear signal on SWI outlining the outer surface of the cortex,\n"
        "within the adjacent subarachnoid space, or both.\n\n"
        "Each row: SWI zoom (40 mm)  |  same with the candidate outlined in red  |\n"
        "8 mm minIP slab (veins appear as continuous branching tubes).\n"
        "Image orientation: neurological (patient left on the left of the image).\n\n"
        "Candidates are in RANDOM order; no scores are shown.\n"
        "Write the codes on the answer grid (last page)."), fontsize=9.5, va="top", family="monospace")
    f.text(0.06, 0.18, "Also: please mark any cSS you see that has NO number (outline it in freeview,\n"
           "or note slice and location) - it is needed to measure what the detector misses.", fontsize=10)
    pdf.savefig(f); plt.close(f)
    nums = sorted(sheet)
    for p0 in range(0, len(nums), PER_PAGE):
        f, axs = plt.subplots(PER_PAGE, 3, figsize=(8.27, 11.69))
        f.subplots_adjust(left=0.08, right=0.98, top=0.96, bottom=0.03, wspace=0.04, hspace=0.12)
        for r in range(PER_PAGE):
            for a in axs[r]: a.axis("off")
            if p0 + r >= len(nums): continue
            n = nums[p0 + r]; m = mask_of(n)
            if not m.any(): continue
            pts = np.argwhere(m); ci, cj = pts[:, 0].mean(), pts[:, 1].mean()
            kc = int(np.clip(np.bincount(pts[:, 2]).argmax(), 0, I.shape[2] - 1))
            i0, i1 = int(max(ci - HW, 0)), int(min(ci + HW, I.shape[0]))
            j0, j1 = int(max(cj - HW, 0)), int(min(cj + HW, I.shape[1]))
            sl = I[i0:i1, j0:j1, kc].T; ms = m[i0:i1, j0:j1, kc].T
            k0, k1 = max(kc - SLAB, 0), min(kc + SLAB + 1, I.shape[2])
            mip = np.where(I[i0:i1, j0:j1, k0:k1] > 0, I[i0:i1, j0:j1, k0:k1], hi).min(axis=2).T
            for c, (img, outline) in enumerate(((sl, False), (sl, True), (mip, True))):
                a = axs[r][c]; a.imshow(img, cmap="gray", origin="lower", vmin=lo, vmax=hi)
                if outline and ms.any():
                    src = m[i0:i1, j0:j1, k0:k1].any(axis=2).T if c == 2 else ms
                    a.contour(src.astype(float), levels=[0.5], colors="red", linewidths=0.9)
            axs[r][0].text(-0.06, 0.5, f"{n}", transform=axs[r][0].transAxes, fontsize=16, weight="bold",
                           ha="right", va="center")
            if r == 0 or p0 + r == p0:
                for c, t in enumerate(("SWI", "outlined", "minIP 8 mm")): axs[r][c].set_title(t, fontsize=9)
        pdf.savefig(f); plt.close(f)
    # answer grid
    f = plt.figure(figsize=(8.27, 11.69)); f.text(0.08, 0.95, f"Answers - {subj}   (codes on page 1, e.g. C, VS:4)", fontsize=14, weight="bold")
    f.text(0.08, 0.925, "Reader: ____________________    Date: ____________", fontsize=10)
    cols = 4; rows = int(np.ceil(len(nums) / cols))
    for n in nums:
        c, r = (n - 1) // rows, (n - 1) % rows
        f.text(0.08 + c * 0.22, 0.89 - r * (0.85 / max(rows, 1)), f"{n:>3}  ______", fontsize=10, family="monospace")
    pdf.savefig(f); plt.close(f)
print(f"{subj}: expert sheet with {len(nums)} objects ({sum(1 for v in sheet.values() if v[0] == 'excluded')} "
      f"excluded mixed in) -> {out}\n"
      f"      key (do NOT give to the expert): {base}/review/{subj}_expert_key.csv")
