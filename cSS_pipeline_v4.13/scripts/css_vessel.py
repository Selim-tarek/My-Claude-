"""Bounded local vein tracking for detect_css.py (v4.16) - RECORDED ONLY, no rule uses these values.

vein_tree_mm saturated on real SWI because the dark tubular network merges brain-wide. These features
stay local and follow the vessel map (work/ID_vessels.nii.gz, the same VES map detect_css.py builds):

  cand_geodesic_mm      length of the candidate along its own voxels (two-sweep geodesic diameter)
  local_cont_max_mm     dark continuation beyond the candidate's ends, the longer of the two ends,
  local_cont_min_mm     ... and the shorter end. Each end is followed geodesically through vessel or dark
                        voxels (dark: so a vein draining into a wide dural sinus continues into it),
                        only forwards (inside a 60 deg cone around the end's outward tangent), up to
                        VES_REACH_MM; 1-voxel gaps are bridged at a cost
  local_cont_total_mm   sum of both ends
  terminal_ends         ends with < 3 mm continuation (2 = stops at both ends, cSS-like; 0 = runs on)
  local_net_mm          skeleton length of the vessel network connected to the candidate within reach
  local_branch_points   junctions in that network with >= 3 arms of >= 3 mm (one count per junction)
  local_branch_per10mm  junctions per 10 mm of network skeleton
  end_type_max/_min     (v4.18) what each end's continuation does, for the end with the longer / shorter
                        continuation: "stops" (< 3 mm), "runs_on" (continues without a junction),
                        "branches" (its continuation contains a junction with >= 3 arms of >= 3 mm),
                        "midline" (reaches the mid-sagittal zone, toward the sagittal sinus / falx, from an
                        end that is not itself there). Two "stops" is the cSS pattern; a vein runs on.
Tracking stays in the pial band (cortex + subarachnoid CSF): the target is veins running ALONG the surface.

A cSS focus touching an unrelated vein inherits that vein's network: read these with the images."""
import numpy as np
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra, connected_components
from skimage.morphology import skeletonize

CONE_DEG = 60.0          # forward cone at each end (veins curve; cSS has no forward continuation at all)
TANGENT_MM = 3.0         # end tangent from the point this far back along the candidate
TERMINAL_MM = 3.0        # continuation below this = the candidate ends there
ARM_MIN_MM = 3.0         # a junction counts when >= 3 arms of at least this length meet
GAP_COST = 1.5           # bridging a 1-voxel gap costs 1.5 x its length

# neighbour offsets up to 2 voxels apart (2 = bridges a 1-voxel gap); forward half only
_OFF = np.array([(a, b, c) for a in range(-2, 3) for b in range(-2, 3) for c in range(-2, 3)
                 if (a, b, c) > (0, 0, 0) and max(abs(a), abs(b), abs(c)) <= 2])


def _graph(idx, shape, vox):
    """sparse graph over voxel coordinates idx (n x 3): 26-neighbours, plus gap-bridging edges"""
    lin = np.ravel_multi_index(idx.T, shape); order = np.argsort(lin); sl = lin[order]
    rows, cols, w = [], [], []
    for o in _OFF:
        q = idx + o
        ok = np.all((q >= 0) & (q < np.array(shape)), axis=1)
        ql = np.full(len(idx), -1, np.int64); ql[ok] = np.ravel_multi_index(q[ok].T, shape)
        pos = np.clip(np.searchsorted(sl, ql), 0, len(sl) - 1)
        hit = ok & (sl[pos] == ql)
        d = float(np.linalg.norm(o * vox)) * (GAP_COST if np.abs(o).max() > 1 else 1.0)
        rows.append(np.nonzero(hit)[0]); cols.append(order[pos[hit]]); w.append(np.full(hit.sum(), d))
    r, c, ww = np.concatenate(rows), np.concatenate(cols), np.concatenate(w)
    n = len(idx)
    return coo_matrix((np.r_[ww, ww], (np.r_[r, c], np.r_[c, r])), shape=(n, n)).tocsr()


def _skeleton_stats(net, vox):
    """skeleton length (mm) and junctions with >= 3 arms of >= ARM_MIN_MM"""
    if net.sum() < 3:
        return 0.0, 0
    sk = skeletonize(net) > 0
    step = float(np.mean(vox))
    length = float(sk.sum()) * step
    k = np.ones((3, 3, 3), int); k[1, 1, 1] = 0
    nb = ndi.convolve(sk.astype(int), k, mode="constant")
    junc = sk & (nb >= 3)
    if not junc.any():
        return length, 0
    jl, nj = ndi.label(junc, structure=np.ones((3, 3, 3), bool))
    arms, na = ndi.label(sk & ~ndi.binary_dilation(junc, structure=np.ones((3, 3, 3), bool)),
                         structure=np.ones((3, 3, 3), bool))
    arm_len = ndi.sum(np.ones_like(arms), arms, index=np.arange(1, na + 1)) * step if na else np.array([])
    long_arm = np.r_[False, arm_len >= ARM_MIN_MM]
    n_valid = 0
    for j, sl in enumerate(ndi.find_objects(jl), 1):
        if sl is None: continue
        big = tuple(slice(max(s.start - 2, 0), s.stop + 2) for s in sl)
        ring = ndi.binary_dilation(jl[big] == j, structure=np.ones((3, 3, 3), bool), iterations=2)
        touching = np.unique(arms[big][ring]); touching = touching[touching > 0]
        n_valid += int(long_arm[touching].sum() >= 3)
    return length, n_valid


def _end_type(cont_mm, reached, end_idx, shape, vox, midline):
    """stops / runs_on / branches / midline for one end (see module docstring)"""
    if cont_mm < TERMINAL_MM or not len(reached):
        return "stops"
    if midline is not None and not midline[tuple(end_idx)] and midline[tuple(reached.T)].any():
        return "midline"
    m = np.zeros(shape, bool); m[tuple(reached.T)] = True; m[tuple(end_idx)] = True
    return "branches" if _skeleton_stats(m, vox)[1] > 0 else "runs_on"


def local_tracking(cm, ves, vox, reach_mm, band=None, dark=None, midline=None):
    """cm: candidate mask, ves: vessel mask (same box, candidate +- reach_mm), vox: voxel size (mm),
    band: where tracking may go (detect_css.py passes the pial band: cortex + CSF, not deep white matter -
    continuation INTO the WM is vessel_wm_mm's job, and WM noise otherwise fakes continuation),
    dark: dark voxels (z < -2); continuation may also run through them, so that a vein draining into a
    dural sinus continues into it (a sinus is too wide to be a Hessian 'tube' at sigma <= 1.6 mm),
    midline: mid-sagittal zone (same box), for end_type midline"""
    path_ok = ves | dark if dark is not None else ves.copy()
    if band is not None:
        ves = ves & band; path_ok &= band
    vox = np.asarray(vox, float)
    nan = dict(cand_geodesic_mm=np.nan, local_cont_max_mm=np.nan, local_cont_min_mm=np.nan,
               local_cont_total_mm=np.nan, terminal_ends=np.nan, local_net_mm=np.nan,
               local_branch_points=np.nan, local_branch_per10mm=np.nan, end_type_max="", end_type_min="")
    ci = np.argwhere(cm)
    if len(ci) < 3:
        return nan
    # ---- candidate geodesic diameter and its two ends (largest connected piece)
    Gc = _graph(ci, cm.shape, vox)
    d0 = dijkstra(Gc, indices=0)
    comp = np.isfinite(d0)
    if comp.sum() < 0.5 * len(ci):                      # start from a node of the biggest piece
        _, lab = connected_components(Gc)
        start = int(np.argmax(lab == np.bincount(lab).argmax())); d0 = dijkstra(Gc, indices=start)
    a = int(np.nanargmax(np.where(np.isfinite(d0), d0, -1)))
    da, pred = dijkstra(Gc, indices=a, return_predecessors=True)
    b = int(np.nanargmax(np.where(np.isfinite(da), da, -1)))
    L = float(da[b]); out = dict(nan, cand_geodesic_mm=round(L, 1))
    path = [b]
    while path[-1] != a and pred[path[-1]] >= 0:
        path.append(int(pred[path[-1]]))
    path = np.array(path[::-1])                         # a ... b along the candidate
    dpath = da[path]
    # ---- forward continuation from each end through vessel voxels only (candidate voxels excluded)
    nodes_v = np.argwhere(path_ok & ~cm)
    cont, etype = [], []
    m_ = min(TANGENT_MM, L / 2)
    back_a = path[min(np.searchsorted(dpath, m_), len(path) - 1)]
    back_b = path[max(np.searchsorted(dpath, L - m_, side="right") - 1, 0)]
    for end, back in ((a, back_a), (b, back_b)):
        pe = ci[end] * vox; t = pe - ci[back] * vox
        if np.linalg.norm(t) < 1e-6 or not len(nodes_v):
            cont.append(0.0); etype.append("stops"); continue
        t /= np.linalg.norm(t)
        rel = nodes_v * vox - pe; dist = np.linalg.norm(rel, axis=1)
        fwd = (dist <= reach_mm) & ((rel @ t >= np.cos(np.radians(CONE_DEG)) * dist) | (dist <= 2.0 * vox.max()))
        if not fwd.any():
            cont.append(0.0); etype.append("stops"); continue
        idx = np.vstack([ci[end][None], nodes_v[fwd]])
        dd = dijkstra(_graph(idx, cm.shape, vox), indices=0, limit=reach_mm)
        c_ = float(np.max(dd[np.isfinite(dd)])); cont.append(c_)
        etype.append(_end_type(c_, idx[1:][np.isfinite(dd[1:])], ci[end], cm.shape, vox, midline))
    order = np.argsort(cont)[::-1]
    out.update(end_type_max=etype[order[0]], end_type_min=etype[order[1]])
    out.update(local_cont_max_mm=round(max(cont), 1), local_cont_min_mm=round(min(cont), 1),
               local_cont_total_mm=round(sum(cont), 1),
               terminal_ends=int(sum(c < TERMINAL_MM for c in cont)))
    # ---- vessel network connected to the candidate within reach: skeleton length and junctions
    allv = np.argwhere(ves | cm)
    if len(allv) >= 3:
        G = _graph(allv, cm.shape, vox)
        src = np.nonzero(cm[tuple(allv.T)])[0]
        dn = dijkstra(G, indices=src, limit=reach_mm, min_only=True)
        net = np.zeros(cm.shape, bool); net[tuple(allv[np.isfinite(dn)].T)] = True
        length, nb = _skeleton_stats(net | cm, vox)
        out.update(local_net_mm=round(length, 1), local_branch_points=int(nb),
                   local_branch_per10mm=round(10.0 * nb / max(length, 1.0), 2))
    return out
