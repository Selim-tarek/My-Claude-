#!/usr/bin/env python3
"""cSS candidate detector v3.
usage: detect_css.py SUBJECT [z_thr=-2.5] [ridge_pct=85] [rim_mm=3]

Changes from v2 (ideas from published tools):
  - SHIVA-CMB: clip bright outliers at the 99.5th percentile before normalising
  - van Harten 2023 / hysteresis: weak voxels (z<-1.5) are kept when connected to a
    strong seed (z<-2.5) -> fewer broken-up lesions, better medium-strength detection
  - microbleednet / Cao 2026: record features that separate cSS from veins
    (cortex fraction, distance from cortex, branching) for review and a future classifier;
    only low cortex fraction (<0.35, typical of mid-sulcus veins) is penalised for now
  - long structures spanning >60% of the slab (sinus / fissure veins) are down-ranked
  - FIX: volumes were under-reported ~2.7x (regionprops already returns mm3 with spacing)
  - MGH cSS histology: 'surface_gradient' = how much darker the pial surface is than
    1-2.5 mm deeper (recorded only; does not change ranking until validated)
  - saves a 'dark' map so scoring can grow accepted lesions to their full extent
"""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
from skimage.filters import sato
from skimage.measure import regionprops
from skimage.morphology import skeletonize

subj   = sys.argv[1]
z_thr  = float(sys.argv[2]) if len(sys.argv) > 2 else -2.5
r_pct  = float(sys.argv[3]) if len(sys.argv) > 3 else 85
RIM_MM = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
Z_LOW, R_PCT_LOW = -2.0, 80           # hysteresis (weak) thresholds (v3.2: tighter, stops leaking into normal cortex)
GROW_Z = -1.5                          # 'dark' map used for region growing at scoring
MIN_MM3, MIN_SLICES, MIN_ELONG = 25.0, 2, 2.5   # 25 mm3 true volume ~= old setting
ARTIFACT_ZONES = ("orbitofrontal", "entorhinal", "temporalpole", "parahippocampal",
                  "fusiform", "inferiortemporal")
LH = [2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 17, 18, 26, 28]
RH = [41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60]

base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
swi = nib.load(f"{base}/data/{subj}_swi.nii")
I   = swi.get_fdata().astype(np.float32)
seg = nib.load(f"{base}/work/{subj}_seg_swispace.nii.gz").get_fdata().astype(np.int32)
vox = tuple(float(v) for v in swi.header.get_zooms()[:3]); vmm3 = float(np.prod(vox))
s2d = np.ones((3, 3, 1), bool); s3 = np.ones((3, 3, 3), bool)

brain  = ndi.binary_fill_holes(I > 0)
cortex = (seg >= 1000) | (seg == 3) | (seg == 42)
it = max(1, int(round(2.0 / vox[0])))
edge  = brain & ~ndi.binary_erosion(brain, structure=s2d, iterations=it)
shell = (ndi.binary_dilation(cortex, structure=s2d, iterations=2) | edge) & brain
surface = (~brain) | (seg == 24)
near_surface = ndi.binary_dilation(surface, structure=s2d, iterations=it)

# 1. robust normalisation
hi = np.percentile(I[brain], 99.5)
Ic = np.minimum(I, hi)
c = Ic[cortex & (I > 0)]
mu = float(np.median(c)); sd = float(1.4826 * np.median(np.abs(c - mu)) + 1e-6)
# v3.2 LOCAL reference: compare each voxel with cortex within ~10 mm, so regionally dark but
# normal cortex (e.g. iron-rich motor cortex) is not called abnormal
_sig = tuple(10.0 / v for v in vox)
_cw = ndi.gaussian_filter(cortex.astype(np.float32), _sig)
_cm = ndi.gaussian_filter((Ic * cortex).astype(np.float32), _sig)
local_mu = np.where(_cw > 0.05, _cm / np.maximum(_cw, 1e-6), mu)
z = (Ic - local_mu) / sd

# 2. skull-strip edge removal (3D, mm)
dist_mm = ndi.distance_transform_edt(brain, sampling=vox)
rim = brain & (dist_mm < RIM_MM)
If = Ic.copy(); If[~brain] = mu

# 3. curvilinear dark structures, slice by slice
R = np.zeros_like(I)
for k in range(I.shape[2]):
    if I[:, :, k].max() > 0:
        R[:, :, k] = sato(If[:, :, k], sigmas=[1, 2, 3], black_ridges=True)
zone = shell & ~rim
r_hi = np.percentile(R[shell], r_pct); r_lo = np.percentile(R[shell], R_PCT_LOW)
strong = zone & (z < z_thr) & (R > r_hi)
weak   = zone & (z < Z_LOW) & (R > r_lo)          # strong is a subset of weak

# 4. hemisphere map; hysteresis labelling separately in each hemisphere
Lm = np.isin(seg, LH) | ((seg >= 1000) & (seg < 2000))
Rm = np.isin(seg, RH) | (seg >= 2000)
hemiL = ndi.distance_transform_edt(~Lm, sampling=vox) <= ndi.distance_transform_edt(~Rm, sampling=vox)
dmid = np.where(hemiL, ndi.distance_transform_edt(hemiL, sampling=vox),
                ndi.distance_transform_edt(~hemiL, sampling=vox))
lab = np.zeros(I.shape, np.int32); n_raw = 0; nxt = 0
for side in (hemiL, ~hemiL):
    wl, nw = ndi.label(weak & side, structure=s3)
    n_raw += nw
    keep = np.unique(wl[strong & side]); keep = keep[keep > 0]
    remap = np.zeros(nw + 1, np.int32)
    remap[keep] = np.arange(1, len(keep) + 1) + nxt
    nz = wl > 0
    lab[nz] = remap[wl[nz]]
    nxt += len(keep)

os.makedirs(f"{base}/work", exist_ok=True)
nib.save(nib.Nifti1Image((zone & (z < GROW_Z)).astype(np.uint8), swi.affine),
         f"{base}/work/{subj}_dark.nii.gz")

# 5. features
cortex_dist = ndi.distance_transform_edt(~cortex, sampling=vox)      # 0 inside cortex
# depth below the pial surface for cortex voxels (distance to nearest non-cortex, non-WM voxel)
pial_depth = ndi.distance_transform_edt(cortex | np.isin(seg, [2, 41]), sampling=vox)
n_brain_sl = max(1, int(brain.any(axis=(0, 1)).sum()))
lut = {}
lutp = os.path.join(os.environ.get("FREESURFER_HOME", ""), "FreeSurferColorLUT.txt")
if os.path.exists(lutp):
    for line in open(lutp):
        t = line.split()
        if len(t) > 2 and t[0].isdigit():
            lut[int(t[0])] = t[1]

def branches_per_10mm(mm):
    """2D skeleton per slice: branch points per 10 mm of centre-line (veins branch, cSS rarely)"""
    bp = 0; length = 0.0
    k8 = np.ones((3, 3)); k8[1, 1] = 0
    for k in range(mm.shape[2]):
        sl = mm[:, :, k]
        if sl.sum() < 3: continue
        sk = skeletonize(sl)
        nb = ndi.convolve(sk.astype(int), k8, mode="constant")
        bp += int((sk & (nb >= 3)).sum()); length += sk.sum() * vox[0]
    return 10.0 * bp / max(length, 1.0)

objs = ndi.find_objects(lab)
rows = []
for p in regionprops(lab, spacing=vox):
    vol = float(p.area)                         # already mm3 because spacing is given
    if vol < MIN_MM3: continue
    sl = objs[p.label - 1]
    pad = tuple(slice(max(s.start - 3, 0), s.stop + 3) for s in sl)
    mm = lab[pad] == p.label
    n_sl = int(mm.any(axis=(0, 1)).sum())
    if n_sl < MIN_SLICES: continue
    elong = p.axis_major_length / max(p.axis_minor_length, 1e-3)
    if elong < MIN_ELONG: continue
    nb = seg[pad][ndi.binary_dilation(mm, iterations=2)]
    nb = nb[nb >= 1000]
    if nb.size == 0: continue
    # hemisphere = majority of the candidate's own voxels in the hemisphere map (v4 fix: v3 counted
    # own L/R labels, but cSS lies in sulcal CSF where both counts are 0, so medial candidates
    # could take the region - and the score - of the opposite hemisphere)
    is_left = bool(hemiL[pad][mm].mean() >= 0.5)
    side = nb[nb < 2000] if is_left else nb[nb >= 2000]
    if side.size: nb = side
    code = int(np.bincount(nb).argmax())
    region = lut.get(code, str(code))
    contact   = float(near_surface[pad][mm].mean())
    depth     = float(-z[pad][mm].mean())
    cfrac     = float(cortex[pad][mm].mean())
    cdist     = float(np.median(cortex_dist[pad][mm]))
    bdens     = branches_per_10mm(mm)
    # depth profile (MGH histology: cSS iron is densest in the outermost cortex):
    # darker at the pial surface than 1-2.5 mm deeper -> positive gradient
    nbh = ndi.binary_dilation(mm, structure=s2d, iterations=4) & cortex[pad]
    dp = pial_depth[pad]; zz = z[pad]
    outer = nbh & (dp <= 1.0); inner = nbh & (dp > 1.0) & (dp <= 2.5)
    sgrad = float(zz[inner].mean() - zz[outer].mean()) if outer.sum() >= 5 and inner.sum() >= 5 else float('nan')
    artz      = any(a in region for a in ARTIFACT_ZONES)
    midz      = bool(np.median(dmid[pad][mm]) < 5.0)
    vein_like = bool(cfrac < 0.35)   # branching is recorded only: thick cSS bands also give spiky skeletons
    longs     = bool(n_sl * vox[2] > 40.0)    # v3.2: > 40 mm top-to-bottom (sinus / large veins)
    score = depth * (0.25 + contact) * np.log1p(vol)
    for flag, f in ((artz, 0.5), (midz, 0.5), (vein_like, 0.6), (longs, 0.5)):
        if flag: score *= f
    c0 = np.argwhere(mm).mean(0) + [s.start for s in pad]
    rows.append(dict(hemi="L" if is_left else "R", region=region, label=code,
                     volume_mm3=round(vol, 1), n_slices=n_sl, elongation=round(elong, 1),
                     surface_contact=round(contact, 2), darkness_z=round(depth, 1),
                     artifact_zone=int(artz), midline_zone=int(midz), score=round(score, 2),
                     vox_i=int(round(c0[0])), vox_j=int(round(c0[1])), vox_k=int(round(c0[2])),
                     cortex_frac=round(cfrac, 2), cortex_dist_mm=round(cdist, 2),
                     branch_per10mm=round(bdens, 2), surface_gradient=round(sgrad, 2),
                     vein_like=int(vein_like),
                     long_structure=int(longs), _lab=p.label))

cols = ["cand_id", "hemi", "region", "label", "volume_mm3", "n_slices", "elongation",
        "surface_contact", "darkness_z", "artifact_zone", "midline_zone", "score",
        "vox_i", "vox_j", "vox_k", "cortex_frac", "cortex_dist_mm", "branch_per10mm", "surface_gradient",
        "vein_like", "long_structure", "accept"]
keep = np.zeros(I.shape, np.int16)
if rows:
    df = pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)
    df.insert(0, "cand_id", np.arange(1, len(df) + 1))
    remap = np.zeros(lab.max() + 1, np.int16)
    remap[df["_lab"].values] = df["cand_id"].values
    keep = remap[lab]
    df["accept"] = ""
    df = df[cols]
else:
    df = pd.DataFrame(columns=cols)

os.makedirs(f"{base}/review", exist_ok=True)
df.to_csv(f"{base}/review/{subj}_candidates.csv", index=False)
nib.save(nib.Nifti1Image(keep, swi.affine), f"{base}/review/{subj}_candidates.nii.gz")
print(f"{subj}: {n_raw} raw components -> {len(df)} candidates "
      f"({int(df.artifact_zone.astype(int).sum()) if len(df) else 0} artifact zone, "
      f"{int(df.midline_zone.astype(int).sum()) if len(df) else 0} midline, "
      f"{int(df.vein_like.astype(int).sum()) if len(df) else 0} vein-like)\n")
if len(df):
    print(df.head(15)[["cand_id", "hemi", "region", "volume_mm3", "n_slices", "darkness_z",
                       "cortex_frac", "branch_per10mm", "vein_like", "score",
                       "vox_i", "vox_j", "vox_k"]].to_string(index=False))
