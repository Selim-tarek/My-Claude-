"""Derive PH6 from the v4 phantom PH3: add the hardest real mimic - a cortical vein running ALONG the pial
surface (0.8 mm outside it, over many slices), with a side branch, draining into a superior sagittal sinus.
usage: phantom_surfvein.py   (needs PH3 in $CSS_BASE; writes PH6 data, labels and work/PH6_surfvein.nii.gz:
       1 = main surface vein, 2 = its branch, 3 = superior sagittal sinus)
PH3's cSS lesions, veins and ICH are unchanged (same truth files are copied)."""
import os, shutil, numpy as np, nibabel as nib
from scipy import ndimage as ndi
B = os.environ["CSS_BASE"]
img = nib.load(f"{B}/data/PH3_swi.nii"); I = img.get_fdata().astype(np.float32)
vox = img.header.get_zooms()[:3]; sh = I.shape
x, y, z = np.meshgrid(*[(np.arange(n) - n / 2) * v for n, v in zip(sh, vox)], indexing="ij")
r = np.sqrt((x / 52) ** 2 + (y / 58) ** 2 + (z / 44) ** 2)
mask = (1 - r) * 52 > -4.5
P = np.stack([x, y, z], -1)

def surf_point(alpha, beta, out_mm):
    """point at polar angle alpha / azimuth beta (deg), out_mm outside the ellipsoid pial surface"""
    a, b = np.radians(alpha), np.radians(beta)
    u = np.array([np.sin(a) * np.cos(b), np.sin(a) * np.sin(b), np.cos(a)])
    s = 1.0 / np.sqrt((u[0] / 52) ** 2 + (u[1] / 58) ** 2 + (u[2] / 44) ** 2)
    return u * s * (1 + out_mm / 52)

def tube(points, radius):
    """voxels within radius of a polyline (sampled every <= 0.3 mm)"""
    m = np.zeros(sh, bool)
    pts = [points[0]]
    for p0, p1 in zip(points[:-1], points[1:]):
        n = max(int(np.ceil(np.linalg.norm(p1 - p0) / 0.3)), 1)
        pts += [p0 + (p1 - p0) * t for t in np.linspace(0, 1, n + 1)[1:]]
    ijk = np.round(np.array(pts) / np.array(vox) + np.array(sh) / 2).astype(int)
    ijk = ijk[np.all((ijk >= 0) & (ijk < np.array(sh)), axis=1)]
    m[tuple(ijk.T)] = True
    return (ndi.distance_transform_edt(~m, sampling=vox) <= radius) & mask

# main vein: from the lateral convexity (polar 75 deg) up to the vertex midline, over the gyral crown
# between sulci 14 and 15 (azimuth 157.5 deg, left hemisphere), 0.8 mm outside the pial surface
azim = lambda a: 157.5 + 4 * np.sin(np.radians(a * 3))      # gentle meander
main_pts = [surf_point(a, azim(a), 0.8) for a in np.linspace(75, 4, 120)]
# side branch leaving the main vein at polar 45 deg (joined to it) and running 20 mm antero-laterally
branch_pts = [surf_point(a, b, 0.8) for a, b in zip(np.linspace(45, 62, 40), np.linspace(azim(45), 138, 40))]
# superior sagittal sinus: midline tube 2.5 mm outside the surface, polar -60..60 deg in the y-z plane
sss_pts = [np.array([0.0, 58 * np.sin(t) * (1 + 2.5 / 52), 44 * np.cos(t) * (1 + 2.5 / 52)])
           for t in np.radians(np.linspace(-60, 60, 200))]
lab = np.zeros(sh, np.int16)
for k, (pts, rad) in enumerate(((main_pts, 0.7), (branch_pts, 0.6), (sss_pts, 1.8)), 1):
    lab[tube(pts, rad) & (lab == 0)] = k
w = np.clip(ndi.gaussian_filter((lab > 0).astype(np.float32), [0.4 / v for v in vox]) * 1.6, 0, 1)
I = I * (1 - 0.6 * w)
I[~mask] = 0; I[mask] = np.maximum(I[mask], 1)
nib.save(nib.Nifti1Image(I, img.affine), f"{B}/data/PH6_swi.nii")
shutil.copy(f"{B}/synthseg/PH3_seg.nii.gz", f"{B}/synthseg/PH6_seg.nii.gz")
for f in ("a2009s_swispace", "truth", "veins", "ich_truth"):
    shutil.copy(f"{B}/work/PH3_{f}.nii.gz", f"{B}/work/PH6_{f}.nii.gz")
nib.save(nib.Nifti1Image(lab, img.affine), f"{B}/work/PH6_surfvein.nii.gz")
print(f"PH6: PH3 + surface-following vein ({int((lab == 1).sum())} vox) with branch ({int((lab == 2).sum())} vox) "
      f"draining to a sagittal sinus ({int((lab == 3).sum())} vox)")
