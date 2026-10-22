#!/usr/bin/env python3
"""Pial-surface (geodesic) features for candidates AND excluded objects - recon-all route only (v4.16).
usage: surface_features.py SUBJECT [--surface-xfm tkras_to_swi.txt] [--near-mm 3] [--reach-mm 25] [--force]
Reads subjects/ID/surf/{lh,rh}.pial, the detector output (review/ID_candidates|excluded.csv/.nii.gz) and
work/ID_dark.nii.gz; writes review/ID_surface_features.csv (one row per object: list, obj_id, fingerprint,
features) and review/ID_surface_features.json (transform, registration check). RECORDED ONLY: no rule,
score or tier uses these values until expert labels show what they are worth.

Surface -> SWI transform. FreeSurfer surfaces are in tkRAS. By default the transform is built from
subjects/ID/mri/orig.mgz (tkRAS -> T1 scanner RAS) and work/ID_swi2t1.lta (prep_anat.sh). The LTA's direction
is NOT trusted: both directions are tried and the one whose pial vertices lie on the label pial boundary of
work/ID_seg_swispace.nii.gz (which mri_vol2vol produced with the same LTA) is used. If neither agrees
within 2 mm the script stops (--force to run anyway). --surface-xfm: a 4x4 text affine tkRAS -> SWI world.

Features (each object is projected onto its nearest pial vertices, <= --near-mm away):
  surf_valid             1 if >= 3 voxels project onto the surface
  surf_normal_dist_mm    median signed distance from the pial surface along the vertex normal
                         (- inside cortex / subpial, + out in the CSF); surf_normal_iqr_mm its spread
  surf_area_mm2          pial area covered (projected vertices + 1 ring)
  surf_pieces            surface-connected pieces of the projection (>= 5 vertices); 2 when it lines both
                         banks of a sulcus (they are close in space but far apart ON the surface)
  surf_band_continuity   share of the projection in its largest surface-connected piece
  geo_length_mm          geodesic length of that piece, measured inside it (two-sweep diameter)
  geo_euclid_ratio       geodesic / straight-line distance between its two ends (>1 = folds with the cortex)
  surf_dark_cont_mm      how far DARK pial surface continues beyond the projection (<= --reach-mm);
                         vertex dark = work/ID_dark.nii.gz within 1.5 mm along its normal
  surf_exit_paths        separate dark paths leaving the projection 3-8 mm out: 0 = isolated band,
                         1 = runs on one way, 2 = runs through, >= 3 = branching
  surf_sulc / surf_curv   (v4.18) median FreeSurfer ?h.sulc (sulcal depth; > 0 = in a sulcus, < 0 = on a gyral
  surf_thickness_mm       crown) and ?h.curv (mean curvature) and ?h.thickness over the projected vertices;
  surf_sulcal_frac        share of projected vertices with sulc > 0. NaN when the file is missing. Both cSS and
                          sulcal veins sit in sulci: context for the ranking model, not a separator on its own.
  (v4.19) DIMENSIONALITY ON THE SURFACE. In one slice a vein and a cSS band both look like a line; on the pial
  surface a vein stays a thin curve (as wide as the vessel) and cSS is a patch (subpial deposit over several mm
  of cortex, across several slices). Measured on the largest surface piece:
  surf_width_mm          projected area (core vertices, no ring) / geodesic length = mean width of the footprint
  surf_cross_mm          2 x the 95th percentile geodesic distance from the footprint's centre line (the
                         diameter path), inside the piece incl. its 1-vertex ring = width at the wide parts
  surf_inradius_mm       largest geodesic distance from the edge of the piece = half the widest part; a curve
                         stays at ~1 ring (mesh spacing ~0.7-1 mm) however long it is, a patch grows
  surf_aspect            surf_cross_mm / geo_length_mm (0 = curve, ~1 = round patch)
  All NaN when the piece is shorter than 2 mm. A vein lying on a cSS patch, or cSS confined to one line of
  cortex, breaks the rule of thumb - recorded only, read with the images.
This replaces the uploaded geodesic_vein_features.py prototype: there, "continuation" was measured on the
bare mesh (always the search limit, 25 mm) and "branches" were mesh vertices with >= 6 neighbours."""
import sys, os, json, argparse, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, FP
from css_surface import apply, surface_to_swi

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("subject"); ap.add_argument("--surface-xfm", default="")
ap.add_argument("--subjects-dir", default="")
ap.add_argument("--near-mm", type=float, default=3.0); ap.add_argument("--reach-mm", type=float, default=25.0)
ap.add_argument("--force", action="store_true")
a = ap.parse_args(); subj = a.subject; base = _base()
sdir = a.subjects_dir or os.environ.get("CSS_SUBJECTS_DIR", f"{base}/subjects")
swi = nib.load(f"{base}/data/{subj}_swi.nii"); aff = swi.affine; inv_aff = np.linalg.inv(aff)
T, tsrc, checks, surfs = surface_to_swi(base, subj, sdir, a.surface_xfm, a.force)
from css_common import update_manifest, sha256
update_manifest(base, subj, "surface_features", dict(transform_source=tsrc, label_agreement_mm=checks,
    pial={h: sha256(f"{sdir}/{subj}/surf/{h}.pial") for h in ("lh", "rh")}, near_mm=a.near_mm, reach_mm=a.reach_mm))

dark = nib.load(f"{base}/work/{subj}_dark.nii.gz").get_fdata() > 0

class Surf:
    def __init__(self, V, F):
        self.V = apply(T, V); self.F = F.astype(np.int64)
        e = np.r_[F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]; e = np.unique(np.sort(e, axis=1), axis=0)
        w = np.linalg.norm(self.V[e[:, 0]] - self.V[e[:, 1]], axis=1); n = len(self.V)
        self.G = coo_matrix((np.r_[w, w], (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(n, n)).tocsr()
        fn = np.cross(self.V[F[:, 1]] - self.V[F[:, 0]], self.V[F[:, 2]] - self.V[F[:, 0]])   # area-weighted
        N = np.zeros_like(self.V)
        for k in range(3): np.add.at(N, F[:, k], fn)
        self.N = N / (np.linalg.norm(N, axis=1, keepdims=True) + 1e-9)
        if np.median(np.einsum("ij,ij->i", self.V - self.V.mean(0), self.N)) < 0: self.N = -self.N   # outward
        self.area = np.zeros(n); np.add.at(self.area, F.ravel(), np.repeat(np.linalg.norm(fn, axis=1) / 6.0, 3))
        self.tree = cKDTree(self.V)
        # vertex dark: work/ID_dark.nii.gz within 1.5 mm along the normal (subpial or just outside)
        dk = np.zeros(n, bool)
        for t in np.arange(-1.5, 1.51, 0.5):
            ijk = np.round(apply(inv_aff, self.V + t * self.N)).astype(int)
            ok = np.all((ijk >= 0) & (ijk < np.array(dark.shape)), axis=1)
            dk[ok] |= dark[tuple(ijk[ok].T)]
        self.dark = dk
S = {"L": Surf(*surfs["lh"]), "R": Surf(*surfs["rh"])}
# v4.18 per-vertex morphometry from recon-all (read only if present)
for h_, k_ in (("lh", "L"), ("rh", "R")):
    S[k_].morph = {}
    for m_ in ("sulc", "curv", "thickness"):
        p_ = f"{sdir}/{subj}/surf/{h_}.{m_}"
        if os.path.exists(p_):
            v_ = nib.freesurfer.read_morph_data(p_)
            if len(v_) == len(S[k_].V): S[k_].morph[m_] = v_

def sub_dist(sf, nodes, src, limit=np.inf):
    """Dijkstra inside the vertex subset `nodes` (bool) from `src` (indices into the full mesh)"""
    idx = np.nonzero(nodes)[0]; pos = -np.ones(len(nodes), int); pos[idx] = np.arange(len(idx))
    d = dijkstra(sf.G[idx][:, idx], indices=pos[src], limit=limit, min_only=True)
    out = np.full(len(nodes), np.inf); out[idx] = d
    return out

def sub_path(sf, nodes, a_, b_):
    """vertex indices of the shortest path a_ -> b_ inside the vertex subset `nodes`"""
    idx = np.nonzero(nodes)[0]; pos = -np.ones(len(nodes), int); pos[idx] = np.arange(len(idx))
    _, pred = dijkstra(sf.G[idx][:, idx], indices=pos[a_], return_predecessors=True)
    path, k = [], pos[b_]
    while k >= 0:
        path.append(idx[k]); k = pred[k]
    return np.array(path[::-1], int)

def features(m, hemi):
    pts = np.argwhere(m); Pw = apply(aff, pts.astype(float))
    sf = S.get(hemi) or min(S.values(), key=lambda s_: np.median(s_.tree.query(Pw)[0]))
    dd, vi = sf.tree.query(Pw)
    ok = dd <= a.near_mm
    r = dict(surf_valid=int(ok.sum() >= 3))
    if not r["surf_valid"]:
        return r
    off = np.einsum("ij,ij->i", Pw[ok] - sf.V[vi[ok]], sf.N[vi[ok]])
    r["surf_normal_dist_mm"] = round(float(np.median(off)), 2)
    r["surf_normal_iqr_mm"] = round(float(np.subtract(*np.percentile(off, [75, 25]))), 2)
    vv = np.unique(vi[ok]); mo = getattr(sf, "morph", {})
    for m_, col in (("sulc", "surf_sulc"), ("curv", "surf_curv"), ("thickness", "surf_thickness_mm")):
        r[col] = round(float(np.median(mo[m_][vv])), 3) if m_ in mo else np.nan
    r["surf_sulcal_frac"] = round(float((mo["sulc"][vv] > 0).mean()), 2) if "sulc" in mo else np.nan
    P = np.zeros(len(sf.V), bool); P[vi[ok]] = True
    P1 = P.copy(); P1[sf.G[np.nonzero(P)[0]].indices] = True          # + 1 ring: voxel grid vs mesh spacing
    r["surf_area_mm2"] = round(float(sf.area[P1].sum()), 1)
    idx = np.nonzero(P1)[0]
    nc, lab = connected_components(sf.G[idx][:, idx], directed=False)
    sizes = np.bincount(lab)
    r["surf_pieces"] = int((sizes >= 5).sum())
    r["surf_band_continuity"] = round(float(sizes.max() / len(idx)), 2)
    big = np.zeros(len(sf.V), bool); big[idx[lab == sizes.argmax()]] = True
    d0 = sub_dist(sf, big, [int(np.nonzero(big)[0][0])])
    e1 = int(np.argmax(np.where(np.isfinite(d0), d0, -1))); d1 = sub_dist(sf, big, [e1])
    e2 = int(np.argmax(np.where(np.isfinite(d1), d1, -1)))
    r["geo_length_mm"] = round(float(d1[e2]), 1)
    eu = float(np.linalg.norm(sf.V[e1] - sf.V[e2]))
    r["geo_euclid_ratio"] = round(float(d1[e2] / eu), 2) if eu > 1.0 else np.nan
    # v4.19 dimensionality on the surface: curve (vein) vs patch (cSS)
    L = float(d1[e2])
    if L >= 2.0:
        r["surf_width_mm"] = round(float(sf.area[P & big].sum() / L), 2)
        dpath = sub_dist(sf, big, sub_path(sf, big, e1, e2))
        r["surf_cross_mm"] = round(float(2 * np.percentile(dpath[big & np.isfinite(dpath)], 95)), 2)
        bi = np.nonzero(big)[0]
        edge = bi[np.array([(~big[sf.G[v].indices]).any() for v in bi])]
        if len(edge) and len(edge) < len(bi):
            de = sub_dist(sf, big, edge)
            r["surf_inradius_mm"] = round(float(de[big & np.isfinite(de)].max()), 2)
        else:
            r["surf_inradius_mm"] = np.nan
        r["surf_aspect"] = round(r["surf_cross_mm"] / L, 2)
    else:
        r.update(surf_width_mm=np.nan, surf_cross_mm=np.nan, surf_inradius_mm=np.nan, surf_aspect=np.nan)
    # dark pial surface beyond the projection, and how many separate dark paths leave it
    D = sf.dark | P1
    dc = sub_dist(sf, D, np.nonzero(P1)[0], limit=a.reach_mm)
    out = np.isfinite(dc) & ~P1
    r["surf_dark_cont_mm"] = round(float(dc[out].max()), 1) if out.any() else 0.0
    ring = out & (dc >= 3.0) & (dc <= 8.0)
    if ring.any():
        ri = np.nonzero(ring)[0]; _, rl = connected_components(sf.G[ri][:, ri], directed=False)
        r["surf_exit_paths"] = int((np.bincount(rl) >= 3).sum())
    else:
        r["surf_exit_paths"] = 0
    return r

rows = []
for lst, idc in (("candidate", "cand_id"), ("excluded", "excl_id")):
    cp, mp = f"{base}/review/{subj}_{lst if lst == 'excluded' else 'candidates'}.csv", \
        f"{base}/review/{subj}_{lst if lst == 'excluded' else 'candidates'}.nii.gz"
    if not os.path.exists(cp): continue
    d = pd.read_csv(cp)
    if not len(d): continue
    M = nib.load(mp).get_fdata().astype(int)
    for _, rr in d.iterrows():
        m = M == int(rr[idc])
        if not m.any(): continue
        rows.append(dict(list=lst, obj_id=int(rr[idc]), **{f: rr[f] for f in FP},
                         **features(m, str(rr.get("hemi", "")))))
cols = ["list", "obj_id"] + FP + ["surf_valid", "surf_normal_dist_mm", "surf_normal_iqr_mm", "surf_area_mm2",
        "surf_pieces", "surf_band_continuity", "geo_length_mm", "geo_euclid_ratio", "surf_dark_cont_mm",
        "surf_exit_paths", "surf_sulc", "surf_curv", "surf_thickness_mm", "surf_sulcal_frac",
        "surf_width_mm", "surf_cross_mm", "surf_inradius_mm", "surf_aspect"]
pd.DataFrame(rows, columns=cols).to_csv(f"{base}/review/{subj}_surface_features.csv", index=False)
json.dump(dict(subject=subj, transform_source=tsrc, label_agreement_mm=checks, tkras_to_swi=T.tolist(),
               near_mm=a.near_mm, reach_mm=a.reach_mm, note="recorded only; not validated"),
          open(f"{base}/review/{subj}_surface_features.json", "w"), indent=2)
print(f"{subj}: surface features for {len(rows)} objects -> {base}/review/{subj}_surface_features.csv")
