"""Derive PH4 from the v4 phantom PH3: the lowest part becomes cerebellum (labels 8/47) and a dark
curvilinear line is drawn along its surface (a normal tentorial / cerebellar-folia line).
usage: phantom_infra.py   (needs PH3 in $CSS_BASE)"""
import os, numpy as np, nibabel as nib
B = os.environ["CSS_BASE"]
img = nib.load(f"{B}/data/PH3_swi.nii"); I = img.get_fdata().astype(np.float32)
seg = nib.load(f"{B}/synthseg/PH3_seg.nii.gz").get_fdata().astype(np.int32)
vox = img.header.get_zooms()[:3]; sh = I.shape
x, y, z = np.meshgrid(*[(np.arange(n) - n / 2) * v for n, v in zip(sh, vox)], indexing="ij")
low = z < -30
tis = (seg > 0) & (seg != 24)
seg[low & tis] = np.where(x < 0, 8, 47)[low & tis]
r = np.sqrt((x / 52) ** 2 + (y / 58) ** 2 + (z / 44) ** 2); depth = (1 - r) * 52
line = low & (z > -38) & (depth > 0) & (depth < 1.2) & (np.abs(np.arctan2(y, x) - np.radians(-100)) < np.radians(12))
I[line] *= 0.35
nib.save(nib.Nifti1Image(I, img.affine), f"{B}/data/PH4_swi.nii")
nib.save(nib.Nifti1Image(seg, img.affine), f"{B}/synthseg/PH4_seg.nii.gz")
nib.save(nib.Nifti1Image(line.astype(np.uint8), img.affine), f"{B}/work/PH4_infra_line.nii.gz")
print(f"PH4: cerebellum below z=-30 mm, dark cerebellar surface line {int(line.sum())} voxels")
