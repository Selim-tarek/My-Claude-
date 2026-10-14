"""Fake FreeSurfer outputs for a phantom so surface_features.py can be tested without FreeSurfer.
usage: phantom_surfaces.py NAME [shift_mm=6] [--lta-vox2vox]
  --lta-vox2vox: write the LTA as FreeSurfer's bbregister does (type 0 LINEAR_VOX_TO_VOX with src/dst volume info
  blocks using the "cras" key) instead of type 1 RAS->RAS
Writes subjects/NAME/surf/{lh,rh}.pial (marching cubes on the phantom's labelled tissue, in tkRAS),
subjects/NAME/mri/orig.mgz (a 'T1' grid whose scanner RAS is the SWI's shifted by shift_mm in x, y and z),
work/NAME_swi2t1.lta (type 1 RAS->RAS, SWI -> T1, i.e. +shift), and work/NAME_tkras_to_swi.txt (the true
tkRAS -> SWI world affine, for the --surface-xfm route), and label/{lh,rh}.aparc.a2009s.annot (the phantom's
Destrieux-like labels sampled at the vertices, names S_kNN / G_kNN)."""
import os, sys, numpy as np, nibabel as nib
from scipy import ndimage as ndi
from skimage.measure import marching_cubes
from nibabel.freesurfer.io import write_geometry, write_annot
B = os.environ["CSS_BASE"]; name = sys.argv[1]
pos = [x for x in sys.argv[2:] if not x.startswith("--")]; shift = float(pos[0]) if pos else 6.0
swi = nib.load(f"{B}/data/{name}_swi.nii"); aff = swi.affine
seg = nib.load(f"{B}/work/{name}_seg_swispace.nii.gz").get_fdata().astype(int)
brain = ndi.binary_fill_holes(np.asarray(swi.dataobj) > 0)
tissue = ndi.binary_fill_holes(brain & (seg > 0) & (seg != 24))
left = np.isin(seg, [2, 3]) | ((seg >= 1000) & (seg < 2000))
hemi_l = ndi.distance_transform_edt(~left) <= ndi.distance_transform_edt(~(tissue & ~left))
R = np.eye(4); R[:3, 3] = shift                                 # SWI scanner RAS -> "T1" scanner RAS
orig_aff = R @ aff
orig = nib.MGHImage(np.zeros(swi.shape, np.uint8), orig_aff)
for d_ in ("surf", "mri", "label"): os.makedirs(f"{B}/subjects/{name}/{d_}", exist_ok=True)
a2 = nib.load(f"{B}/work/{name}_a2009s_swispace.nii.gz").get_fdata().astype(int)
_, near_a2 = ndi.distance_transform_edt(a2 == 0, return_indices=True)        # nearest labelled cortex voxel
nib.save(orig, f"{B}/subjects/{name}/mri/orig.mgz")
tkr2t1 = orig.header.get_vox2ras() @ np.linalg.inv(orig.header.get_vox2ras_tkr())
swi_to_tkr = np.linalg.inv(tkr2t1) @ R
for h, hm in (("lh", hemi_l), ("rh", ~hemi_l)):
    v, f, _, _ = marching_cubes(np.pad((tissue & hm).astype(np.float32), 1), 0.5)
    v -= 1                                                        # undo the padding
    vw = v @ aff[:3, :3].T + aff[:3, 3]                            # SWI scanner RAS
    vt = vw @ swi_to_tkr[:3, :3].T + swi_to_tkr[:3, 3]             # tkRAS of the fake T1
    write_geometry(f"{B}/subjects/{name}/surf/{h}.pial", vt, f.astype(np.int32))
    # v4.18 morphometry with known values: sulc = height below the top of the hemisphere (mm) - its median,
    # curv = 0.1, thickness = 2.5 mm (check_surfvein.py recomputes the expected medians)
    sulc = (vw[:, 2].max() - vw[:, 2]); sulc = sulc - np.median(sulc)
    nib.freesurfer.write_morph_data(f"{B}/subjects/{name}/surf/{h}.sulc", sulc.astype(np.float32))
    nib.freesurfer.write_morph_data(f"{B}/subjects/{name}/surf/{h}.curv", np.full(len(vt), 0.1, np.float32))
    nib.freesurfer.write_morph_data(f"{B}/subjects/{name}/surf/{h}.thickness", np.full(len(vt), 2.5, np.float32))
    ijk = np.clip(np.round(v).astype(int), 0, np.array(a2.shape) - 1)
    codes = a2[tuple(n_[tuple(ijk.T)] for n_ in near_a2)]
    uc = np.unique(codes)
    names = [("S_k%02d" % (c % 100 - 50) if c % 100 >= 50 else "G_k%02d" % (c % 100 - 1)) if c else "Unknown"
             for c in uc]
    ctab = np.array([[(37 * k) % 255 + 1, (91 * k) % 255, (53 * k) % 255, 0, 0] for k in range(len(uc))], np.int32)
    ctab[:, 4] = ctab[:, 0] + ctab[:, 1] * 256 + ctab[:, 2] * 65536
    write_annot(f"{B}/subjects/{name}/label/{h}.aparc.a2009s.annot", np.searchsorted(uc, codes), ctab, names,
                fill_ctab=False)
def volinfo(tag, A, shape, fn):
    vs = np.linalg.norm(A[:3, :3], axis=0); D = A[:3, :3] / vs
    c = A[:3, :3] @ (np.array(shape[:3]) / 2.0) + A[:3, 3]
    return (f"{tag}\nvalid = 1  # volume info valid\nfilename = {fn}\nvolume = {shape[0]} {shape[1]} {shape[2]}\n"
            f"voxelsize = {vs[0]:.15e} {vs[1]:.15e} {vs[2]:.15e}\n"
            + "".join(f"{k}   = {D[0, i]:.15e} {D[1, i]:.15e} {D[2, i]:.15e}\n" for i, k in enumerate(("xras", "yras", "zras")))
            + f"cras   = {c[0]:.15e} {c[1]:.15e} {c[2]:.15e}\n")
with open(f"{B}/work/{name}_swi2t1.lta", "w") as fh:
    if "--lta-vox2vox" in sys.argv:
        V = np.linalg.inv(orig_aff) @ R @ aff                     # SWI voxel -> T1 voxel
        fh.write("# transform file (phantom, bbregister-style)\ntype      = 0 # LINEAR_VOX_TO_VOX\nnxforms   = 1\n"
                 "mean      = 0.0000 0.0000 0.0000\nsigma     = 1.0000\n1 4 4\n"
                 + "\n".join(" ".join(f"{x:.15e}" for x in row) for row in V) + "\n"
                 + volinfo("src volume info", aff, swi.shape, f"{B}/data/{name}_swi.nii")
                 + volinfo("dst volume info", orig_aff, swi.shape, f"{B}/subjects/{name}/mri/orig.mgz")
                 + "subject phantom\nfscale 0.100000\n")
    else:
        fh.write("# transform file (phantom)\ntype      = 1 # LINEAR_RAS_TO_RAS\nnxforms   = 1\nmean      = 0.0 0.0 0.0\n"
                 "sigma     = 1.0\n1 4 4\n" + "\n".join(" ".join(f"{x:.6f}" for x in row) for row in R) + "\n")
np.savetxt(f"{B}/work/{name}_tkras_to_swi.txt", np.linalg.inv(swi_to_tkr))
print(f"{name}: phantom pial surfaces + orig.mgz + LTA (shift {shift} mm) in {B}/subjects/{name}")
