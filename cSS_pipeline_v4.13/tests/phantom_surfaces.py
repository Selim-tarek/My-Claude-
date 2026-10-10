"""Fake FreeSurfer outputs for a phantom so surface_features.py can be tested without FreeSurfer.
usage: phantom_surfaces.py NAME [shift_mm=6]
Writes subjects/NAME/surf/{lh,rh}.pial (marching cubes on the phantom's labelled tissue, in tkRAS),
subjects/NAME/mri/orig.mgz (a 'T1' grid whose scanner RAS is the SWI's shifted by shift_mm in x, y and z),
work/NAME_swi2t1.lta (type 1 RAS->RAS, SWI -> T1, i.e. +shift), and work/NAME_tkras_to_swi.txt (the true
tkRAS -> SWI world affine, for the --surface-xfm route)."""
import os, sys, numpy as np, nibabel as nib
from scipy import ndimage as ndi
from skimage.measure import marching_cubes
from nibabel.freesurfer.io import write_geometry
B = os.environ["CSS_BASE"]; name = sys.argv[1]; shift = float(sys.argv[2]) if len(sys.argv) > 2 else 6.0
swi = nib.load(f"{B}/data/{name}_swi.nii"); aff = swi.affine
seg = nib.load(f"{B}/work/{name}_seg_swispace.nii.gz").get_fdata().astype(int)
brain = ndi.binary_fill_holes(np.asarray(swi.dataobj) > 0)
tissue = ndi.binary_fill_holes(brain & (seg > 0) & (seg != 24))
left = np.isin(seg, [2, 3]) | ((seg >= 1000) & (seg < 2000))
hemi_l = ndi.distance_transform_edt(~left) <= ndi.distance_transform_edt(~(tissue & ~left))
R = np.eye(4); R[:3, 3] = shift                                 # SWI scanner RAS -> "T1" scanner RAS
orig_aff = R @ aff
orig = nib.MGHImage(np.zeros(swi.shape, np.uint8), orig_aff)
os.makedirs(f"{B}/subjects/{name}/surf", exist_ok=True); os.makedirs(f"{B}/subjects/{name}/mri", exist_ok=True)
nib.save(orig, f"{B}/subjects/{name}/mri/orig.mgz")
tkr2t1 = orig.header.get_vox2ras() @ np.linalg.inv(orig.header.get_vox2ras_tkr())
swi_to_tkr = np.linalg.inv(tkr2t1) @ R
for h, hm in (("lh", hemi_l), ("rh", ~hemi_l)):
    v, f, _, _ = marching_cubes(np.pad((tissue & hm).astype(np.float32), 1), 0.5)
    v -= 1                                                        # undo the padding
    vw = v @ aff[:3, :3].T + aff[:3, 3]                            # SWI scanner RAS
    vt = vw @ swi_to_tkr[:3, :3].T + swi_to_tkr[:3, 3]             # tkRAS of the fake T1
    write_geometry(f"{B}/subjects/{name}/surf/{h}.pial", vt, f.astype(np.int32))
with open(f"{B}/work/{name}_swi2t1.lta", "w") as fh:
    fh.write("# transform file (phantom)\ntype      = 1 # LINEAR_RAS_TO_RAS\nnxforms   = 1\nmean      = 0.0 0.0 0.0\n"
             "sigma     = 1.0\n1 4 4\n" + "\n".join(" ".join(f"{x:.6f}" for x in row) for row in R) + "\n")
np.savetxt(f"{B}/work/{name}_tkras_to_swi.txt", np.linalg.inv(swi_to_tkr))
print(f"{name}: phantom pial surfaces + orig.mgz + LTA (shift {shift} mm) in {B}/subjects/{name}")
