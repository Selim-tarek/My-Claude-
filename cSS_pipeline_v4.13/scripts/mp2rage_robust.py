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
is closer to the measured |INV1|. --lambda auto (v4.18.2): the smallest of 1, 2, 4 ... 128 that leaves < 5 % of the
background band 5-20 mm outside the head brighter than 25 % of the brain, without changing the brain (r >= 0.95). Output is scaled to 0..4095 like a scanner UNI. Prints the
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
# v4.18.2: noise from the background band 5-20 mm OUTSIDE the head (INV2 has a clean head outline). The image border
# is often zero-padded / defaced (LEMON: border noise 10 -> lambda 1, speckle left around the skull and neck).
from skimage.filters import threshold_otsu
vox = np.array(iu.header.get_zooms()[:3], float)
pos = A2[A2 > 0]
head = A2 > threshold_otsu(pos) * 0.5
head = ndi.binary_closing(head, iterations=2); head = ndi.binary_fill_holes(head)
lab, n = ndi.label(head)
if n: head = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
d_out = ndi.distance_transform_edt(~head, sampling=vox)
shell = (d_out > 5) & (d_out <= 20) & (A2 > 0)
if shell.sum() < 1000:                                    # tiny field of view: fall back to everything outside
    shell = (d_out > 2) & (A2 > 0)
noise = float(np.mean(A2[shell])) if shell.any() else float(np.percentile(pos, 1))
core = ndi.binary_erosion(head, iterations=max(1, int(round(10 / vox.min()))))
if not core.any(): core = head
ref = (u[core] + 0.5) * 4095
def robust(lam):
    b = lam * noise ** 2
    return (I1 * A2 - b) / (I1 ** 2 + A2 ** 2 + 2 * b)
lams = [1, 2, 4, 8, 16, 32, 64, 128] if a.lam == "auto" else [float(a.lam)]
prev = None
for lam in lams:
    out_l = np.clip((robust(lam) + 0.5) * 4095, 0, 4095)
    hm = float(np.median(out_l[core]))
    speck = float((out_l[shell] > 0.25 * hm).mean()) if shell.any() else 0.0
    rc = float(np.corrcoef(out_l[core], ref)[0, 1])
    if rc < 0.95 and prev is not None:                    # brain contrast would change: keep the previous lambda
        lam, out, speck, ok = prev; break
    out, ok = out_l, speck < 0.05
    prev = (lam, out_l, speck, ok)
    if ok: break
bd = shell
c = float(np.corrcoef(out[core], ref)[0, 1])
img = nib.Nifti1Image(out.astype(np.float32), iu.affine, iu.header); img.set_data_dtype(np.float32)
nib.save(img, a.out)
print(f"{a.out}: robust MP2RAGE combination, lambda {lam:g} (background noise {noise:.1f}, band 5-20 mm outside the head); "
      f"speckle left in that band {speck:.1%} (was {float(((u[shell] + 0.5) * 4095 > 0.25 * np.median(ref)).mean()):.1%}); "
      f"correlation with the original UNI inside the brain {c:.3f}" + ("" if ok else "  WARNING: background still speckled"))
if c < 0.95: print("WARNING: brain contrast changed (correlation < 0.95) - check the image before recon-all")
