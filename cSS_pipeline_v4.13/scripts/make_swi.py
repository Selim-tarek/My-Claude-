#!/usr/bin/env python3
"""Build an SWI image from gradient-echo MAGNITUDE + PHASE (v4.18) - for datasets that store the raw GRE
(e.g. OpenNeuro ds000221 / MPI-Leipzig LEMON: acq-mag_GRE + acq-phase_GRE) instead of a processed SWI.
usage: make_swi.py MAG.nii.gz PHASE.nii.gz OUT_swi.nii.gz [--sigma-mm 4] [--power 4] [--phase-sign auto|+1|-1]
                   [--hp-out OUT_hpphase.nii.gz]
Standard SWI processing (Haacke 2004): phase -> radians, homodyne high-pass per axial slice
(angle(z * conj(lowpass z)), Gaussian low-pass sigma --sigma-mm, wrap-safe; the same filter the reviewer
uses), a phase mask that darkens PARAMAGNETIC voxels (veins, haemosiderin) linearly to 0 at pi, and
SWI = magnitude x mask^power. Nothing else is changed: the output keeps the magnitude's grid and affine.

Phase sign: which sign is paramagnetic differs between vendors and reconstructions. With --phase-sign auto
(default) it is read from the image: paramagnetic voxels are also dark on the magnitude at this TE, so the
sign of the median high-pass phase of the darkest-for-their-neighbourhood brain voxels is used. The
evidence is printed; if it is weak the script stops and asks for --phase-sign. Check the result once per
protocol: veins must be DARKER on the SWI than on the magnitude (calcium is not darkened).
Research software, not for clinical reading."""
import sys, os, argparse, numpy as np, nibabel as nib
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import phase_to_radians, homodyne_highpass

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("mag"); ap.add_argument("phase"); ap.add_argument("out")
ap.add_argument("--sigma-mm", type=float, default=4.0); ap.add_argument("--power", type=int, default=4)
ap.add_argument("--phase-sign", default="auto"); ap.add_argument("--hp-out", default="")
a = ap.parse_args()

mi, pi_ = nib.load(a.mag), nib.load(a.phase)
if mi.shape[:3] != pi_.shape[:3] or not np.allclose(mi.affine, pi_.affine, atol=1e-3):
    sys.exit("magnitude and phase are not on the same grid - they must come from the same acquisition")
M = np.asarray(mi.dataobj, np.float32); P = np.asarray(pi_.dataobj, np.float32)
if M.ndim == 4 or P.ndim == 4:
    sys.exit("4-D (multi-echo) input: pass one echo (e.g. the longest usable TE) as 3-D images")
vox = np.array(mi.header.get_zooms()[:3], float)
si = int(np.argmax(np.abs(mi.affine[2, :3])))          # superior-inferior array axis: slice-wise filter
P, unit = phase_to_radians(P)                            # scanner units (e.g. -4096..4095) -> radians
# homodyne high-pass, slice by slice in the plane perpendicular to the S-I axis (v5.0: shared with the detector)
HP = homodyne_highpass(P, vox, si, a.sigma_mm)

# brain-ish mask for the sign estimate: bright magnitude, eroded away from the head edge
thr = 0.15 * float(np.percentile(M[M > 0], 99)) if (M > 0).any() else 0.0
head = ndi.binary_fill_holes(M > thr)
lab, n = ndi.label(head)
if n: head = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
core = ndi.binary_erosion(head, iterations=max(1, int(round(10 / vox.min()))))
if a.phase_sign == "auto":
    rel = M / (ndi.gaussian_filter(M, 3.0 / vox) + 1e-6)
    dark = core & (rel < 0.75) & (np.abs(HP) > 0.2)
    if dark.sum() < 200:
        sys.exit(f"too few dark brain voxels ({int(dark.sum())}) to read the phase sign - give --phase-sign +1 or -1")
    med = float(np.median(HP[dark])); frac = float((np.sign(HP[dark]) == np.sign(med)).mean())
    if abs(med) < 0.05 or frac < 0.6:
        sys.exit(f"phase sign unclear (median {med:+.2f} rad, {frac:.0%} agree) - give --phase-sign +1 or -1")
    sign = 1 if med > 0 else -1
    how = f"auto: dark-on-magnitude voxels median high-pass phase {med:+.2f} rad, {frac:.0%} same sign, n={int(dark.sum())}"
else:
    sign = int(a.phase_sign); how = "given"
    if sign not in (1, -1): sys.exit("--phase-sign must be auto, +1 or -1")

# phase mask: paramagnetic side (sign) goes linearly from 1 at 0 to 0 at pi
para = sign * HP
mask = np.where(para > 0, np.clip(1.0 - para / np.pi, 0.0, 1.0), 1.0).astype(np.float32)
S = (M * mask ** a.power).astype(np.float32)
out = nib.Nifti1Image(S, mi.affine, mi.header); out.set_data_dtype(np.float32)
nib.save(out, a.out)
if a.hp_out:
    nib.save(nib.Nifti1Image(HP, mi.affine), a.hp_out)
red = float(1 - np.median(S[core & (mask < 0.9)] / (M[core & (mask < 0.9)] + 1e-6))) if (core & (mask < 0.9)).any() else 0.0
print(f"{os.path.basename(a.out)}: SWI = magnitude x phase-mask^{a.power}; phase {unit}; high-pass sigma {a.sigma_mm} mm "
      f"in-plane (S-I axis {si}); paramagnetic sign {sign:+d} ({how}); {(core & (mask < 0.9)).mean() if core.any() else 0:.1%} "
      f"of brain voxels darkened, median by {red:.0%}")
