#!/usr/bin/env python3
"""Insert synthetic cSS into a (presumed cSS-negative) host scan -> {HOST}S.
usage: make_synthetic.py HOST [n_lesions=8] [depth=0.6] [seed=1]
v2 (ideas from Cao et al. 2026): per-lesion random darkness, thickness (1-2 voxels),
length 10-30 mm over 2-4 slices, random edge blur, and added noise inside lesions,
so lesions are less 'clean' than before."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib
from scipy import ndimage as ndi

host  = sys.argv[1]
N     = int(sys.argv[2])   if len(sys.argv) > 2 else 8
depth = float(sys.argv[3]) if len(sys.argv) > 3 else 0.6
seed  = int(sys.argv[4])   if len(sys.argv) > 4 else 1
rng = np.random.default_rng(seed)
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))

img = nib.load(f"{base}/data/{host}_swi.nii"); I = img.get_fdata().astype(np.float32)
seg = nib.load(f"{base}/work/{host}_seg_swispace.nii.gz").get_fdata().astype(int)
vox = np.array(img.header.get_zooms()[:3], float)
s2 = np.ones((3, 3, 1), bool)
brain = ndi.binary_fill_holes(I > 0)
dist = ndi.distance_transform_edt(brain, sampling=vox)
cortex = seg >= 1000
CONVEX = [3, 8, 11, 15, 22, 24, 27, 28, 29, 30, 31]
convex = cortex & np.isin(seg % 1000, CONVEX)
sulcal = brain & ~cortex & ~np.isin(seg, [2, 41])
boundary = convex & ndi.binary_dilation(sulcal, structure=s2) & (dist > 4.0)
c = I[cortex & (I > 0)]; sd_c = 1.4826 * np.median(np.abs(c - np.median(c)))

nz = I.shape[2]
pts_all = np.argwhere(boundary)
pts_all = pts_all[(pts_all[:, 2] >= 3) & (pts_all[:, 2] < nz - 5)]
if len(pts_all) == 0: sys.exit("no admissible cortical boundary found")

truth, seeds, tries = np.zeros(I.shape, np.int16), [], 0
dark = np.zeros(I.shape, np.float32)
while len(seeds) < N and tries < 8000:
    tries += 1
    p = pts_all[rng.integers(len(pts_all))]
    if any(np.linalg.norm((p - q) * vox) < 25 for q in seeds): continue
    L_mm, nsl = rng.uniform(10, 30), int(rng.integers(2, 5))
    lesion, cur, ok = np.zeros(I.shape, bool), p[:2], True
    for dk in range(nsl):
        k = p[2] + dk
        if k >= nz: ok = False; break
        b = boundary[:, :, k]; cand = np.argwhere(b)
        if len(cand) == 0: ok = False; break
        d = np.linalg.norm((cand - cur) * vox[:2], axis=1); j = d.argmin()
        if d[j] > 4: ok = False; break
        start = tuple(cand[j]); cur = cand[j]
        target, seen, queue = int(L_mm / vox[0]), {start}, [start]
        while queue and len(seen) < target:
            x, y = queue.pop(0)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    nx_, ny_ = x + dx, y + dy
                    if 0 <= nx_ < b.shape[0] and 0 <= ny_ < b.shape[1] and b[nx_, ny_] and (nx_, ny_) not in seen:
                        seen.add((nx_, ny_)); queue.append((nx_, ny_))
        if len(seen) < 0.6 * target: ok = False; break
        for x, y in seen: lesion[x, y, k] = True
    if not ok: continue
    lesion = ndi.binary_dilation(lesion, structure=s2, iterations=int(rng.integers(1, 3))) & brain
    if (truth[lesion] > 0).any(): continue
    seeds.append(p); truth[lesion] = len(seeds)
    d_i = float(np.clip(depth + rng.uniform(-0.1, 0.1), 0.05, 0.9))
    w = ndi.gaussian_filter(lesion.astype(np.float32), sigma=(rng.uniform(0.6, 1.2),) * 2 + (0,))
    dark = np.maximum(dark, np.clip(w * 1.5, 0, 1) * d_i)

Is = I * (1 - dark)
Is += (rng.normal(0, 0.35 * sd_c, I.shape) * (dark > 0.05)).astype(np.float32)
Is[brain] = np.maximum(Is[brain], 1.0)        # keep the brain mask unchanged
Is[~brain] = I[~brain]

out = nib.Nifti1Image(Is.astype(np.float32), img.affine, img.header); out.set_data_dtype(np.float32)
nib.save(out, f"{base}/data/{host}S_swi.nii")
nib.save(nib.Nifti1Image(truth, img.affine), f"{base}/work/{host}S_truth.nii.gz")

lut = {}
lutp = os.path.join(os.environ.get("FREESURFER_HOME", ""), "FreeSurferColorLUT.txt")
if os.path.exists(lutp):
    for line in open(lutp):
        t = line.split()
        if len(t) > 2 and t[0].isdigit(): lut[int(t[0])] = t[1]
print(f"{host}S: inserted {len(seeds)} synthetic cSS lesions (depth ~{depth}, seed {seed})")
for i in range(1, len(seeds) + 1):
    m = truth == i; s = seg[m][seg[m] >= 1000]
    lab = int(np.bincount(s).argmax()) if s.size else 0
    ks = np.unique(np.argwhere(m)[:, 2]); c0 = np.argwhere(m).mean(0).round().astype(int)
    print(f"  lesion {i}: {lut.get(lab, lab)!s:30s} slices {ks.min()}-{ks.max()}  vox {c0[0]} {c0[1]} {c0[2]}")
