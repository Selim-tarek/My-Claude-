"""FreeSurfer pial surface -> SWI space, with a registration check (v4.16; shared by surface_features.py and
surface_overlay.py). FreeSurfer surfaces are in tkRAS. The transform is built from subjects/ID/mri/orig.mgz
(tkRAS -> T1 scanner RAS) and work/ID_swi2t1.lta (prep_anat.sh). The LTA's direction is NOT trusted: both
directions are tried and the one whose pial vertices lie on the label pial boundary of
work/ID_seg_swispace.nii.gz (which mri_vol2vol made with the same LTA) is used. surface_to_swi() stops the
script when even the better direction is more than AGREE_MM off, unless force=True."""
import os, sys, numpy as np, nibabel as nib
from scipy import ndimage as ndi
from nibabel.freesurfer.io import read_geometry

AGREE_MM = 2.0


def apply(T, X):
    return X @ T[:3, :3].T + T[:3, 3]


def read_lta(p):
    """4x4 RAS->RAS matrix of an LTA (type 1 RAS2RAS, or type 0 VOX2VOX converted with its volume info)"""
    L = open(p).read().splitlines()
    typ = int(next(l for l in L if l.strip().startswith("type")).split("=")[1].split()[0])
    i = next(k for k, l in enumerate(L) if l.strip().split()[:3] == ["1", "4", "4"])
    M = np.array([[float(x) for x in L[i + 1 + r].split()] for r in range(4)])
    if typ == 1:
        return M
    def v2r(tag):                                     # vox2ras from "src/dst volume info" block
        k = next(j for j, l in enumerate(L) if l.strip().startswith(tag))
        info = {}
        for l in L[k + 1:k + 10]:
            if "=" in l:
                key, val = l.split("=", 1); info[key.strip()] = val.split()
        dims = np.array([float(x) for x in info["volume"][:3]]); vs = np.array([float(x) for x in info["voxelsize"]])
        R = np.c_[[float(x) for x in info["xras"]], [float(x) for x in info["yras"]], [float(x) for x in info["zras"]]] * vs
        c = np.array([float(x) for x in info["c_ras"]])
        A = np.eye(4); A[:3, :3] = R; A[:3, 3] = c - R @ (dims / 2.0)
        return A
    return v2r("dst volume info") @ M @ np.linalg.inv(v2r("src volume info"))


def surface_to_swi(base, subj, sdir, surface_xfm="", force=False):
    """-> (T tkRAS->SWI world, source description, {option: median vertex-to-label-boundary mm}, surfs)
    surfs = {"lh": (vertices_tkRAS, faces), "rh": ...} of the pial surfaces"""
    swi = nib.load(f"{base}/data/{subj}_swi.nii"); vox = np.array(swi.header.get_zooms()[:3], float)
    seg = nib.load(f"{base}/work/{subj}_seg_swispace.nii.gz").get_fdata().astype(np.int32)
    brain = ndi.binary_fill_holes(np.asarray(swi.dataobj) > 0)
    tissue = ndi.binary_fill_holes(brain & (seg > 0) & (seg != 24))
    pial_sd = ndi.distance_transform_edt(~tissue, sampling=vox) - ndi.distance_transform_edt(tissue, sampling=vox)
    inv_aff = np.linalg.inv(swi.affine)
    agree = lambda Vw: float(np.median(np.abs(ndi.map_coordinates(pial_sd, apply(inv_aff, Vw).T, order=1, mode="nearest"))))
    surfs = {h: read_geometry(f"{sdir}/{subj}/surf/{h}.pial") for h in ("lh", "rh")}
    Vall = np.vstack([surfs["lh"][0], surfs["rh"][0]])[::50]
    if surface_xfm:
        T = np.loadtxt(surface_xfm).reshape(4, 4); src = "explicit --surface-xfm"
        checks = {"explicit": agree(apply(T, Vall))}
    else:
        orig = nib.load(f"{sdir}/{subj}/mri/orig.mgz")
        tkr2t1 = orig.header.get_vox2ras() @ np.linalg.inv(orig.header.get_vox2ras_tkr())
        R = read_lta(f"{base}/work/{subj}_swi2t1.lta")
        opts = {"LTA maps SWI->T1": np.linalg.inv(R) @ tkr2t1, "LTA maps T1->SWI": R @ tkr2t1}
        checks = {k: agree(apply(Tk, Vall)) for k, Tk in opts.items()}
        src = min(checks, key=checks.get); T = opts[src]
    print(f"{subj}: surface -> SWI ({src}); median distance of pial vertices to the label pial boundary: "
          + ", ".join(f"{k} {v:.2f} mm" for k, v in checks.items()))
    if min(checks.values()) > AGREE_MM and not force:
        sys.exit(f"{subj}: surfaces do not match the SWI labels (best {min(checks.values()):.2f} mm > {AGREE_MM} mm)"
                 f" - check the registration / transform (or --force)")
    return T, src, checks, surfs
