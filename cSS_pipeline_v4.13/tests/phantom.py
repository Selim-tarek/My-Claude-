"""Synthetic phantom brain (ellipsoid, radial sulci, midline fissure, dark veins) for tests.
usage: phantom.py NAME SEED   (writes into $CSS_BASE)"""
import os, numpy as np, nibabel as nib, sys
B = os.environ["CSS_BASE"]
from scipy import ndimage as ndi
name, sd = sys.argv[1], int(sys.argv[2]); rng = np.random.default_rng(sd)
sh = (220, 220, 40); vox = (0.5, 0.5, 2.0)
aff = np.diag(list(vox) + [1.0])
x, y, zz = np.meshgrid(*[(np.arange(n) - n/2) * v for n, v in zip(sh, vox)], indexing="ij")
r = np.sqrt((x/50)**2 + (y/52)**2 + (zz/38)**2)          # ellipsoid "brain"
brain = r < 1
seg = np.zeros(sh, np.int32)
depth_mm = (1 - r) * 50                                    # approx depth from surface
left = x < 0; front = y > 0
ctxcode = np.where(left, 1000, 2000) + np.where(front, 28, 22)
seg[brain] = np.where(left, 2, 41)[brain]
seg[brain & (depth_mm < 3)] = ctxcode[brain & (depth_mm < 3)]
ang = np.arctan2(y, x)
I = np.zeros(sh, np.float32)
for a0 in np.arange(-np.pi, np.pi, np.pi/9):              # radial sulci
    d_line = np.abs(np.sin(ang - a0)) * np.sqrt(x**2 + y**2)
    sul = brain & (d_line < 0.8) & (depth_mm < 12) & (np.cos(ang - a0) > 0)
    bank = brain & (d_line < 3.8) & (depth_mm < 15) & (np.cos(ang - a0) > 0) & ~sul
    seg[bank] = ctxcode[bank]; seg[sul] = 24
seg[np.abs(x) < 0.7] = np.where(brain[np.abs(x) < 0.7], 24, 0)   # midline fissure
I[np.isin(seg, [2, 41])] = 230; I[seg >= 1000] = 200; I[seg == 24] = 255
I[brain] += rng.normal(0, 12, brain.sum())
# veins: dark tubes along some sulci, with a branch
for a0 in np.arange(-np.pi, np.pi, np.pi/4):
    d_line = np.abs(np.sin(ang - a0 + 0.17)) * np.sqrt(x**2 + y**2)
    v = brain & (d_line < 0.6) & (depth_mm > 2) & (depth_mm < 9) & (np.cos(ang - a0) > 0) & (np.abs(zz) < 14)
    I[v] = 70
I[~brain] = 0; I[brain] = np.maximum(I[brain], 1)
seg[~brain] = 0
nib.save(nib.Nifti1Image(I, aff), f"{B}/data/{name}_swi.nii")
nib.save(nib.Nifti1Image(seg, aff), f"{B}/synthseg/{name}_seg.nii.gz")
