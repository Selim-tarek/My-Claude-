#!/usr/bin/env python3
"""MP2RAGE background clean-up before FreeSurfer (v4.18): the 'robust combination' of O'Brien et al. 2014
(PLoS One 9:e99676) as implemented by J. Marques (RobustCombination.m), from MAGNITUDE images.
usage: mp2rage_robust.py UNI.nii.gz INV1.nii.gz INV2.nii.gz OUT_T1.nii.gz [--lambda auto|1|5|...]
The MP2RAGE UNI image has salt-and-pepper noise outside the head (division of two noisy signals). recon-all
then fails to strip dura / background (LEMON P201: 'extra large defect ... dura has not been removed from
wm.mgz', topology fixing aborted). The robust combination adds a noise term beta to the UNI formula:
      UNI_r = (INV1 * INV2 - beta) / (INV1^2 + INV2^2 + 2 beta),   beta = lambda * mean(INV2 noise)^2
which leaves the brain (high INV2 signal) unchanged and pushes the background to -0.5 (black). INV1's sign
(lost in magnitude images) is recovered from UNI: INV1 = the root of UNI = INV1 INV2 / (INV1^2 + INV2^2) that
is closer to the measured |INV1|. --lambda auto: the smallest of 1, 2, 4, 8, 16, 32 that makes the image
border darker than 10 % of the head's median. Output is scaled to 0..4095 like a scanner UNI. Prints the
correlation with the original UNI inside the head (should stay > 0.95)."""
import sys, argparse, numpy as np, nibabel as nib
from scipy import ndimage as ndi

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("uni"); ap.add_argument("inv1"); ap.add_argument("inv2"); ap.add_argument("out")
ap.add_argument("--lambda", dest="lam", default="auto")
a = ap.parse_args()
iu, i1, i2 = (nib.load(p) for p in (a.uni, a.inv1, a.inv2))
if not (iu.shape[:3] == i1.shape[:3] == i2.shape[:3]) or not (np.allclose(iu.affine, i1.affine, atol=1e-2)
                                                             and np.allclose(iu.affine, i2.affine, atol=1e-2)):
    sys.exit("UNI, INV1 and INV2 must be on the same grid (same MP2RAGE acquisition)")
U = np.asarray(iu.dataobj, np.float64); A1 = np.abs(np.asarray(i1.dataobj, np.float64)); A2 = np.abs(np.asarray(i2.dataobj, np.float64))
# UNI to [-0.5, 0.5]
if U.max() > 1.5:
    lo, hi = (0.0, 4095.0) if U.min() >= 0 and U.max() <= 4096 else (float(U.min()), float(U.max()))
    u = (U - lo) / (hi - lo) - 0.5
else:
    u = U - 0.5 if U.min() >= 0 else U
u = np.clip(u, -0.5, 0.5)
# INV1 sign from UNI: u (i1^2 + i2^2) = i1 i2  ->  u i1^2 - i2 i1 + u i2^2 = 0
with np.errstate(invalid="ignore", divide="ignore"):
    disc = np.sqrt(np.clip(A2 ** 2 - 4 * u * u * A2 ** 2, 0, None))
    r_pos = np.where(np.abs(u) > 1e-6, (A2 + disc) / (2 * u), 0.0)
    r_neg = np.where(np.abs(u) > 1e-6, (A2 - disc) / (2 * u), 0.0)
I1 = np.where(np.abs(np.abs(r_pos) - A1) < np.abs(np.abs(r_neg) - A1), r_pos, r_neg)
I1 = np.where(np.isfinite(I1), I1, 0.0)
# noise of INV2 from the image border (defaced / zero-filled voxels ignored)
bd = np.zeros(A2.shape, bool); w = max(3, min(A2.shape) // 20)
for ax in range(3):
    sl = [slice(None)] * 3; sl[ax] = slice(0, w); bd[tuple(sl)] = True; sl[ax] = slice(-w, None); bd[tuple(sl)] = True
noise = A2[bd & (A2 > 0)]
noise = float(np.mean(noise)) if noise.size else float(np.percentile(A2[A2 > 0], 1))
head = ndi.binary_fill_holes(ndi.binary_opening(A2 > 5 * noise, iterations=2))
lab, n = ndi.label(head)
if n: head = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
def robust(lam):
    b = lam * noise ** 2
    return (I1 * A2 - b) / (I1 ** 2 + A2 ** 2 + 2 * b)
lams = [1, 2, 4, 8, 16, 32] if a.lam == "auto" else [float(a.lam)]
for lam in lams:
    r = robust(lam); out = np.clip((r + 0.5) * 4095, 0, 4095)
    ok = np.median(out[bd & ~head]) < 0.1 * np.median(out[head]) if (bd & ~head).any() else True
    if ok: break
c = float(np.corrcoef(out[head], (u[head] + 0.5) * 4095)[0, 1])
img = nib.Nifti1Image(out.astype(np.float32), iu.affine, iu.header); img.set_data_dtype(np.float32)
nib.save(img, a.out)
print(f"{a.out}: robust MP2RAGE combination, lambda {lam:g} (noise {noise:.1f}); border median "
      f"{np.median(out[bd & ~head]) if (bd & ~head).any() else 0:.0f} vs head {np.median(out[head]):.0f}; "
      f"correlation with the original UNI inside the head {c:.3f}" + ("" if ok else "  WARNING: background still bright"))
if c < 0.95: print("WARNING: brain contrast changed (correlation < 0.95) - check the image before recon-all")
