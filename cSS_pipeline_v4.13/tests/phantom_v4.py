"""Phantom brain with the structures the cSS definition separates (v4 tests).
usage: phantom_v4.py NAME SEED        (writes into $CSS_BASE; appends test labels to the LUT)

Geometry (mm, 0.8 x 0.8 x 1.2 voxels like a clinical SWI):
  ellipsoid hemispheres, 2.5 mm cortex, 4.5 mm extra-axial CSF (SynthStrip keeps CSF), midline
  fissure, 16 radial sulci (2 mm CSF, 14 mm deep).
Truth  work/NAME_truth.nii.gz  (cSS, ids 1-6):
  1-3 tram-track: subpial band on BOTH banks of a sulcus        (the classic sign)
  4-5 convexity: subpial band following the gyral crown          (curvilinear sheet)
  6   single bank of a sulcus
Mimics work/NAME_veins.nii.gz (ids 1-6):
  1-4 tubular veins in the middle of the sulcal CSF (1,2 and 3,4 are left/right mirror pairs)
  5-6 cortical surface veins running tangentially in the CSF just outside the pial surface
  7   transcortical / medullary vein: radial tube from the surface CSF through the cortex 14 mm into WM
Also: a lobar ICH (dark sphere, 8 mm) in left white matter -> work/NAME_ich_truth.nii.gz
Labels: synthseg/NAME_seg.nii.gz (aseg/DK-like) and work/NAME_a2009s_swispace.nii.gz
(Destrieux-like: S_k = banks of sulcus k, G_k = gyrus between sulci k and k+1)."""
import os, sys, numpy as np, nibabel as nib
from scipy import ndimage as ndi

B = os.environ["CSS_BASE"]
name, sd = sys.argv[1], int(sys.argv[2]); rng = np.random.default_rng(sd)
vox = (0.8, 0.8, 1.2); sh = (150, 165, 90)
aff = np.diag(list(vox) + [1.0])
x, y, z = np.meshgrid(*[(np.arange(n) - n / 2) * v for n, v in zip(sh, vox)], indexing="ij")
R0 = 52.0
r = np.sqrt((x / 52) ** 2 + (y / 58) ** 2 + (z / 44) ** 2)
depth = (1 - r) * R0                        # ~mm below the pial surface (negative = outside)
inside = r < 1                              # pial surface
mask = depth > -4.5                         # brain mask incl. 4.5 mm extra-axial CSF (atrophic CAA)
left = x < 0
rho = np.sqrt(x ** 2 + y ** 2); ang = np.arctan2(y, x)

seg = np.zeros(sh, np.int32)
seg[mask] = 24
seg[inside] = np.where(left, 2, 41)[inside]
ctx = np.where(left, 1000, 2000) + np.where(y > 0, 28, 22)
cort = inside & (depth < 2.5)
A0 = np.arange(-np.pi, np.pi, np.pi / 8) + np.pi / 16          # 16 sulci, none on the midline
sul_any = np.zeros(sh, bool); bank_k = np.full(sh, -1, np.int16)
dl = np.zeros((len(A0),) + sh, np.float32)
for k, a0 in enumerate(A0):
    d_line = np.abs(np.sin(ang - a0)) * rho
    front = np.cos(ang - a0) > 0
    dl[k] = np.where(front, d_line, 99)
    sul = inside & (dl[k] < 1.0) & (depth < 14)
    bank = inside & (dl[k] < 3.5) & (depth < 16.5) & ~sul
    cort |= bank; sul_any |= sul
    bank_k[bank & (bank_k < 0)] = k
fiss = mask & (np.abs(x) < 1.0)
cort |= inside & (np.abs(x) < 3.5) & ~fiss
cort &= ~sul_any & ~fiss
seg[cort] = ctx[cort]; seg[sul_any | fiss] = 24
seg[~mask] = 0

# Destrieux-like labels on the cortex
sector = (np.floor((ang - A0[0]) / (np.pi / 8)).astype(int)) % len(A0)
a2 = np.zeros(sh, np.int32)
hemi_off = np.where(left, 11100, 12100)
a2[cort] = (hemi_off + 1 + sector)[cort]                         # G_k
bk = cort & (bank_k >= 0)
a2[bk] = (hemi_off + 50 + bank_k)[bk]                            # S_k

I = np.zeros(sh, np.float32)
I[np.isin(seg, [2, 41])] = 180; I[seg >= 1000] = 210; I[seg == 24] = 250
I[mask] += rng.normal(0, 10, int(mask.sum()))

def darken(m, f=0.6, blur=0.6):
    """hemosiderin / deoxy-Hb: multiplicative darkening with blooming (gaussian edge)"""
    w = np.clip(ndi.gaussian_filter(m.astype(np.float32), [blur / v for v in vox]) * 1.6, 0, 1)
    I[:] = I * (1 - f * w)

truth = np.zeros(sh, np.int16); veins = np.zeros(sh, np.int16)
zc = lambda lo, hi: (z > lo) & (z < hi)
# 1-3 tram-track: subpial 1.2 mm band on both banks (sulci 2, 6 left side; 11 right side)
for i, (k, z0) in enumerate([(5, -6), (7, 4), (13, -10)], 1):
    m = (dl[k] >= 1.0) & (dl[k] < 2.2) & (depth > 3) & (depth < 12) & zc(z0, z0 + 16) & cort
    truth[m] = i
# 4-5 convexity crown: outer 1 mm of cortex over a 14 deg arc, 18 mm high
for i, (a, z0) in enumerate([(A0[4] + np.pi / 16, -2), (A0[12] + np.pi / 16, 6)], 4):
    m = cort & (depth >= 0) & (depth < 1.2) & (np.abs(np.angle(np.exp(1j * (ang - a)))) < np.radians(7)) & zc(z0, z0 + 18)
    truth[m] = i
# 6 single bank
k = 9; s = np.sign(np.sin(ang - A0[k]))
m = (dl[k] >= 1.0) & (dl[k] < 2.2) & (s > 0) & (depth > 3) & (depth < 12) & zc(-4, 12) & cort
truth[m] = 6
# veins 1-4: tubes (r 0.7 mm) in the sulcal CSF at 5 mm depth, running 30 mm along z;
# 1/2 and 3/4 are left-right mirror images (normal anatomy is roughly symmetric)
# (sulcus k mirrors sulcus 7-k mod 16; sulci 5, 7, 9, 13 carry cSS)
for i, k in enumerate([1, 6, 11, 12], 1):
    c_xy = np.array([np.cos(A0[k]), np.sin(A0[k])])
    # point on the sulcal plane where depth ~5 mm (search along the ray)
    t = np.linspace(20, 60, 400); px, py = c_xy[0] * t, c_xy[1] * t
    rr = np.sqrt((px / 52) ** 2 + (py / 58) ** 2); p = np.argmin(np.abs((1 - rr) * R0 - 5))
    m = (np.sqrt((x - px[p]) ** 2 + (y - py[p]) ** 2) < 0.7) & zc(-14, 16) & mask
    veins[m] = i
# veins 5-6: cortical surface veins 1.2 mm outside the pial surface, tangential, 40 deg arc
for i, (a, z0) in enumerate([(np.radians(60), 8), (np.radians(-120), -6)], 5):
    m = (np.abs(depth + 1.2) < 0.7) & (np.abs(z - z0) < 0.7) & \
        (np.abs(np.angle(np.exp(1j * (ang - a)))) < np.radians(20)) & mask
    veins[m] = i
# vein 7: radial (perpendicular to the cortex), gyral crown between sulci 1 and 2, from 2 mm outside
# the pial surface to 14 mm deep; slight tilt so it is seen on several slices
a7 = A0[1] + np.pi / 16
u7 = np.array([np.cos(a7) * 52, np.sin(a7) * 58, 6.0]); u7 /= np.linalg.norm(u7)
P = np.stack([x, y, z - 2.0], -1); tp = P @ u7
dperp = np.linalg.norm(P - tp[..., None] * u7, axis=-1)
m = (dperp < 1.0) & (depth > -2) & (depth < 14) & mask
veins[m] = 7
ich = (np.sqrt((x + 22) ** 2 + (y - 10) ** 2 + (z - 4) ** 2) < 8) & inside
darken(truth > 0, 0.6); darken(veins > 0, 0.6, blur=0.4); darken(ich, 0.7, blur=0.8)
I[~mask] = 0; I[mask] = np.maximum(I[mask], 1)

os.makedirs(f"{B}/data", exist_ok=True); os.makedirs(f"{B}/work", exist_ok=True)
os.makedirs(f"{B}/synthseg", exist_ok=True)
nib.save(nib.Nifti1Image(I, aff), f"{B}/data/{name}_swi.nii")
nib.save(nib.Nifti1Image(seg, aff), f"{B}/synthseg/{name}_seg.nii.gz")
nib.save(nib.Nifti1Image(a2, aff), f"{B}/work/{name}_a2009s_swispace.nii.gz")
nib.save(nib.Nifti1Image(truth, aff), f"{B}/work/{name}_truth.nii.gz")
nib.save(nib.Nifti1Image(veins, aff), f"{B}/work/{name}_veins.nii.gz")
nib.save(nib.Nifti1Image(ich.astype(np.uint8), aff), f"{B}/work/{name}_ich_truth.nii.gz")
lutp = os.path.join(os.environ.get("FREESURFER_HOME", ""), "FreeSurferColorLUT.txt")
if os.path.exists(lutp) and "ctx_lh_S_k00" not in open(lutp).read():
    with open(lutp, "a") as f:
        for h, off in (("lh", 11100), ("rh", 12100)):
            for k in range(len(A0)):
                f.write(f"{off + 1 + k} ctx_{h}_G_k{k:02d} 100 100 100 0\n")
                f.write(f"{off + 50 + k} ctx_{h}_S_k{k:02d} 200 200 200 0\n")
print(f"{name}: phantom v4  cSS lesions {int(truth.max())}  veins {int(veins.max())}  "
      f"ICH {int(ich.sum() * np.prod(vox))} mm3")
