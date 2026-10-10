#!/usr/bin/env python3
"""Insert synthetic cSS (and optional synthetic vein decoys) into a presumed cSS-negative host -> {HOST}S.
usage: make_synthetic.py HOST [n_lesions=8] [depth=0.6] [seed=1] [--tram 0.4] [--veins 4] [--radial 2] [--legacy]

v3 (realistic lesions; after the P006 stress test showed v2 lesions were too big and blob-like, so
volume / branching / elongation separated them for the wrong reason):
  - thin CURVILINEAR lesions: a 1-voxel path traced along the cortex-CSF boundary (no flood fill),
    optionally thickened by 1 voxel (blooming), 5-30 mm long (log-uniform), on 2-5 consecutive slices
  - tram-track: with probability --tram the opposite bank of the same sulcus gets a matching line
  - darkness per lesion around `depth`, edge blur, and noise inside the lesion
  - --veins K: K tubular decoys (radius 0.5-0.9 mm, 15-40 mm) in the sulcal CSF, as dark as the
    lesions; written to work/HOSTS_veins.nii.gz (Cao et al. 2026: vessel decoys reduce false positives)
  - --radial K (v4.11): K transcortical / medullary vein decoys: straight tubes starting 2 mm out in the
    CSF and running perpendicular to the cortex 10-18 mm into the white matter (same veins file)
--legacy reproduces the v2 generator (flood-filled, thicker lesions) for comparison."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib
from scipy import ndimage as ndi

flag = lambda f, d: type(d)(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d
TRAM, NVEIN, NRAD, LEGACY = flag("--tram", 0.4), flag("--veins", 0), flag("--radial", 0), "--legacy" in sys.argv
# positional values given after a flag (e.g. "--veins 4") must not be read as lesion arguments
pos = [a for i, a in enumerate(sys.argv[1:], 1) if not a.startswith("--") and not sys.argv[i - 1] in ("--tram", "--veins", "--radial")]
host = pos[0]; N = int(pos[1]) if len(pos) > 1 else 8
depth = float(pos[2]) if len(pos) > 2 else 0.6; seed = int(pos[3]) if len(pos) > 3 else 1
rng = np.random.default_rng(seed)
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))

img = nib.load(f"{base}/data/{host}_swi.nii"); I = img.get_fdata().astype(np.float32)
seg = nib.load(f"{base}/work/{host}_seg_swispace.nii.gz").get_fdata().astype(int)
vox = np.array(img.header.get_zooms()[:3], float)
s2 = np.ones((3, 3, 1), bool)
brain = ndi.binary_fill_holes(I > 0)
dist = ndi.distance_transform_edt(brain, sampling=vox)
cortex = seg >= 1000
CONVEX = [3, 8, 11, 15, 22, 24, 27, 28, 29, 30, 31]
convex = cortex & np.isin(seg % 1000, CONVEX)
sulcal = brain & ~cortex & ~np.isin(seg, [2, 41])
boundary = convex & ndi.binary_dilation(sulcal, structure=s2) & (dist > 4.0)
c = I[cortex & (I > 0)]; sd_c = 1.4826 * np.median(np.abs(c - np.median(c)))

nz = I.shape[2]
pts_all = np.argwhere(boundary)
pts_all = pts_all[(pts_all[:, 2] >= 3) & (pts_all[:, 2] < nz - 6)]
if len(pts_all) == 0: sys.exit("no admissible cortical boundary found")
NB8 = [(dx, dy) for dx in (-1, 0, 1) for dy in (-1, 0, 1) if dx or dy]

def trace(b, start, L_mm, d0=None):
    """1-voxel path along the boundary b (2-D), preferring to keep its direction (no branches)"""
    path, seen, cur, d, length = [start], {start}, start, d0, 0.0
    while length < L_mm:
        nxt = [(cur[0] + dx, cur[1] + dy) for dx, dy in NB8
               if 0 <= cur[0] + dx < b.shape[0] and 0 <= cur[1] + dy < b.shape[1]
               and b[cur[0] + dx, cur[1] + dy] and (cur[0] + dx, cur[1] + dy) not in seen]
        if not nxt: break
        if d is None: n = nxt[rng.integers(len(nxt))]
        else: n = max(nxt, key=lambda q: (q[0] - cur[0]) * d[0] + (q[1] - cur[1]) * d[1])
        d = (n[0] - cur[0], n[1] - cur[1])
        length += float(np.hypot(d[0] * vox[0], d[1] * vox[1]))
        path.append(n); seen.add(n); cur = n
    return path, length, d

def opposite_bank(path, k):
    """boundary voxels across the sulcal CSF from the path (1.5-6 mm away, CSF in between)"""
    b = boundary[:, :, k]; P = np.array(path)
    cand = np.argwhere(b)
    if not len(cand): return None
    dd = np.sqrt((((cand[:, None, :] - P[None, :, :]) * vox[:2]) ** 2).sum(-1)).min(1)
    cand = cand[(dd > 1.5) & (dd < 6.0)]
    ok = [q for q in map(tuple, cand)
          if sulcal[(q[0] + P[np.argmin(((P - q) ** 2).sum(1))][0]) // 2,
                    (q[1] + P[np.argmin(((P - q) ** 2).sum(1))][1]) // 2, k]]
    return ok[rng.integers(len(ok))] if ok else None

truth, seeds, tries = np.zeros(I.shape, np.int16), [], 0
dark = np.zeros(I.shape, np.float32)
while len(seeds) < N and tries < 8000:
    tries += 1
    p = pts_all[rng.integers(len(pts_all))]
    if any(np.linalg.norm((p - q) * vox) < 25 for q in seeds): continue
    lesion, ok = np.zeros(I.shape, bool), True
    if LEGACY:                                   # v2: flood fill, then 1-2 dilations
        L_mm, nsl, cur = rng.uniform(10, 30), int(rng.integers(2, 5)), p[:2]
        for dk in range(nsl):
            k = p[2] + dk; b = boundary[:, :, k]; cand = np.argwhere(b)
            if k >= nz or len(cand) == 0: ok = False; break
            d = np.linalg.norm((cand - cur) * vox[:2], axis=1); j = d.argmin()
            if d[j] > 4: ok = False; break
            start = tuple(cand[j]); cur = cand[j]
            target, seen, queue = int(L_mm / vox[0]), {start}, [start]
            while queue and len(seen) < target:
                x, y = queue.pop(0)
                for dx, dy in NB8:
                    q = (x + dx, y + dy)
                    if 0 <= q[0] < b.shape[0] and 0 <= q[1] < b.shape[1] and b[q] and q not in seen:
                        seen.add(q); queue.append(q)
            if len(seen) < 0.6 * target: ok = False; break
            for q in seen: lesion[q[0], q[1], k] = True
        if not ok: continue
        lesion = ndi.binary_dilation(lesion, structure=s2, iterations=int(rng.integers(1, 3))) & brain
    else:                                        # v3: thin traced line(s)
        L_mm = float(np.exp(rng.uniform(np.log(5), np.log(30))))
        nsl = int(rng.integers(2, 6)); cur, d0 = tuple(p[:2]), None
        tram = rng.random() < TRAM; opp = None
        for dk in range(nsl):
            k = p[2] + dk; b = boundary[:, :, k]; cand = np.argwhere(b)
            if k >= nz or not len(cand): ok = False; break
            dd = np.linalg.norm((cand - np.array(cur)) * vox[:2], axis=1); j = dd.argmin()
            if dd[j] > 3: ok = False; break
            path, length, d0 = trace(b, tuple(cand[j]), L_mm, d0 if dk else None)
            if length < 0.6 * L_mm: ok = False; break
            for q in path: lesion[q[0], q[1], k] = True
            cur = tuple(cand[j])
            if tram:
                o = opposite_bank(path, k) if opp is None else None
                start = o if o is not None else opp
                if start is not None:
                    bo = boundary[:, :, k]; cand_o = np.argwhere(bo)
                    jo = np.linalg.norm((cand_o - np.array(start)) * vox[:2], axis=1).argmin()
                    po, _, _ = trace(bo, tuple(cand_o[jo]), L_mm)
                    for q in po: lesion[q[0], q[1], k] = True
                    opp = tuple(cand_o[jo])
        if not ok: continue
        if rng.random() < 0.5:                   # blooming: 1-voxel in-plane thickening, stays at the surface
            lesion = ndi.binary_dilation(lesion, structure=s2) & (cortex | sulcal) & brain
    if (truth[lesion] > 0).any() or not lesion.any(): continue
    seeds.append(p); truth[lesion] = len(seeds)
    d_i = float(np.clip(depth + rng.uniform(-0.1, 0.1), 0.05, 0.9))
    w = ndi.gaussian_filter(lesion.astype(np.float32), sigma=(rng.uniform(0.4, 0.8),) * 2 + (0,))
    dark = np.maximum(dark, np.clip(w * 2.0, 0, 1) * d_i)

# ---- vein decoys: tubes running in the sulcal CSF
veins = np.zeros(I.shape, np.int16)
if NVEIN and not LEGACY:
    mid = sulcal & (ndi.distance_transform_edt(sulcal, sampling=vox) >= 0.7) & (dist > 6.0)
    in_csf = ndi.binary_dilation(sulcal, structure=np.ones((3, 3, 3), bool))   # decoys stay in sulcal CSF
    vp = np.argwhere(mid); nv = 0; vt = 0
    while nv < NVEIN and vt < 2000 and len(vp):
        vt += 1
        q = vp[rng.integers(len(vp))].astype(float)
        if (truth[tuple(q.astype(int))] > 0) or any(np.linalg.norm((q - s) * vox) < 15 for s in seeds): continue
        dvec = rng.normal(size=3); dvec[2] *= rng.uniform(0.2, 1.0); dvec /= np.linalg.norm(dvec * vox) + 1e-9
        L, r, pts, pnt = rng.uniform(15, 40), rng.uniform(0.5, 0.9), [], q.copy()
        for _ in range(int(L / 0.5)):
            dvec += rng.normal(scale=0.05, size=3) / vox; dvec /= np.linalg.norm(dvec * vox) + 1e-9
            pnt = pnt + dvec * 0.5
            ii = np.round(pnt).astype(int)
            if (ii < 0).any() or (ii >= np.array(I.shape)).any() or not in_csf[tuple(ii)]: break
            pts.append(ii)
        if len(pts) * 0.5 < 10: continue
        m = np.zeros(I.shape, bool); m[tuple(np.array(pts).T)] = True
        m = (ndi.distance_transform_edt(~m, sampling=vox) <= r) & brain
        if (truth[m] > 0).any() or (veins[m] > 0).any(): continue
        nv += 1; veins[m] = nv
        dark = np.maximum(dark, m * float(np.clip(depth + rng.uniform(-0.1, 0.1), 0.05, 0.9)))

# ---- transcortical / medullary vein decoys: straight tubes perpendicular to the cortex into the WM
if NRAD and not LEGACY:
    wm = np.isin(seg, [2, 41])
    g = np.gradient(ndi.gaussian_filter(ndi.distance_transform_edt(~wm, sampling=vox), 1.0), *vox)
    nv, vt = int(veins.max()), 0; n0 = nv
    while nv - n0 < NRAD and vt < 2000:
        vt += 1
        q = pts_all[rng.integers(len(pts_all))].astype(float)
        if any(np.linalg.norm((q - s_) * vox) < 15 for s_ in seeds) or veins[tuple(q.astype(int))]: continue
        u = -np.array([gg[tuple(q.astype(int))] for gg in g]); nu = np.linalg.norm(u)
        if nu < 0.3: continue
        u = u / nu                                         # mm direction toward the white matter
        L, r = rng.uniform(10, 18), rng.uniform(0.5, 0.9)
        tt = np.arange(-2.0, 3.0 + L, 0.3)                 # 2 mm outside, ~3 mm cortex, L mm into WM
        ii = np.round(q[None, :] + tt[:, None] * u[None, :] / vox).astype(int)
        if (ii < 0).any() or (ii >= np.array(I.shape)).any(): continue
        if wm[tuple(ii[tt > 4.0].T)].mean() < 0.6 or not brain[tuple(ii.T)].all(): continue
        m = np.zeros(I.shape, bool); m[tuple(ii.T)] = True
        m = (ndi.distance_transform_edt(~m, sampling=vox) <= r) & brain
        if (truth[m] > 0).any() or (veins[m] > 0).any(): continue
        nv += 1; veins[m] = nv
        dark = np.maximum(dark, m * float(np.clip(depth + rng.uniform(-0.1, 0.1), 0.05, 0.9)))

Is = I * (1 - dark)
Is += (rng.normal(0, 0.35 * sd_c, I.shape) * (dark > 0.05)).astype(np.float32)
Is[brain] = np.maximum(Is[brain], 1.0)        # keep the brain mask unchanged
Is[~brain] = I[~brain]

out = nib.Nifti1Image(Is.astype(np.float32), img.affine, img.header); out.set_data_dtype(np.float32)
nib.save(out, f"{base}/data/{host}S_swi.nii")
nib.save(nib.Nifti1Image(truth, img.affine), f"{base}/work/{host}S_truth.nii.gz")
nib.save(nib.Nifti1Image(veins, img.affine), f"{base}/work/{host}S_veins.nii.gz")

lut = {}
lutp = os.path.join(os.environ.get("FREESURFER_HOME", ""), "FreeSurferColorLUT.txt")
if os.path.exists(lutp):
    for line in open(lutp):
        t = line.split()
        if len(t) > 2 and t[0].isdigit(): lut[int(t[0])] = t[1]
print(f"{host}S: inserted {len(seeds)} synthetic cSS lesions ({'legacy' if LEGACY else 'v3 thin'}, depth ~{depth}, "
      f"seed {seed}){f', {int(veins.max())} vein decoys' if veins.max() else ''}")
for i in range(1, len(seeds) + 1):
    m = truth == i; s = seg[m][seg[m] >= 1000]
    lab = int(np.bincount(s).argmax()) if s.size else 0
    ks = np.unique(np.argwhere(m)[:, 2]); c0 = np.argwhere(m).mean(0).round().astype(int)
    print(f"  lesion {i}: {lut.get(lab, lab)!s:30s} slices {ks.min()}-{ks.max()}  vox {c0[0]} {c0[1]} {c0[2]}"
          f"  {m.sum() * np.prod(vox):.0f} mm3")
