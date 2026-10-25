"""Feature families of docs/FEATURE_SPEC (cSS vs VS / VC / VT / CA / MB) that the detector did not yet compute (v5.1.2).
Pure functions on a candidate crop; imported by detect_css.py, unit-tested in tests/check_v512.py. RECORDED ONLY: no
rule, score or tier uses them. All distances in mm (voxel size passed explicitly).

Already computed elsewhere under another name (not duplicated): pial_dist_mm (median signed pial distance),
cortex_frac (cortex overlap), parenchyma_frac, bank_frac / tram_frac, tube_ratio, cand_radius_med/p95_mm / cv,
local_cont_max/min_mm (curved continuation), vessel_run_mm (straight continuation), vessel_wm_mm (white-matter
continuation), local_branch_per10mm, terminal_ends, end_type_* ("midline" = reaches the sagittal sinus zone), surf_*
(surface morphology, surface_features.py), phase_z / phase_r (dipole fit), mb_cluster_n, flair_csf_z / flair_ctx_z,
ich_dist_mm. Not computable here: CT hyperdensity (no CT input)."""
import numpy as np
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

FEATURE_SCHEMA = "css-features/1.0"


def csf_maps(csf, vox):
    """whole-brain maps for position_features: distance into the CSF from tissue (0 in tissue), and 'within 1.5 mm
    of CSF'. Computed once per scan (a per-candidate crop would cut the distances at the crop edge)."""
    return (ndi.distance_transform_edt(csf, sampling=vox).astype(np.float32),
            ndi.distance_transform_edt(~csf, sampling=vox) <= 1.5)


def position_features(mm, pial_sd, csf_depth, near_csf, sulcal_depth=None):
    """Family B/D: pial-distance percentiles, contact within 1 / 2 mm of the pial boundary, CSF adjacency, how far
    into the CSF the candidate sits (median distance from tissue: ~0 on a bank, larger in mid-sulcus = VC) and how
    deep below the gyral crowns (sulcal_depth_map). All maps are crops at the candidate."""
    psd = pial_sd[mm]
    out = dict(pial_sd_p10_mm=round(float(np.percentile(psd, 10)), 2), pial_sd_p90_mm=round(float(np.percentile(psd, 90)), 2),
               surface_contact_frac_1mm=round(float((np.abs(psd) <= 1.0).mean()), 2),
               surface_contact_frac_2mm=round(float((np.abs(psd) <= 2.0).mean()), 2),
               csf_adjacent_frac=round(float(near_csf[mm].mean()), 2),
               csf_depth_mm=round(float(np.median(csf_depth[mm])), 2),
               sulcal_depth_mm=round(float(np.median(sulcal_depth[mm])), 2) if sulcal_depth is not None else np.nan)
    return out


def compactness_features(mm, vox):
    """Family G: equivalent-sphere diameter, 3-D compactness 36 pi V^2 / S^3 (S = anisotropic marching-cubes area;
    1 = sphere), sphericity = compactness^(1/3), elongation (longest / shortest principal axis) and flatness
    (middle / shortest), principal axes from the voxel coordinates in mm."""
    vox = np.asarray(vox, float); V = float(mm.sum() * np.prod(vox))
    out = dict(equivalent_sphere_diameter_mm=round(float((6 * V / np.pi) ** (1 / 3)), 2),
               compactness_3d=np.nan, sphericity=np.nan, elongation_3d=np.nan, flatness_3d=np.nan)
    if mm.sum() < 4: return out
    try:
        from skimage.measure import marching_cubes, mesh_surface_area
        v, f, _, _ = marching_cubes(np.pad(mm, 1).astype(np.float32), 0.5, spacing=tuple(vox))
        S = float(mesh_surface_area(v, f))
        if S > 0:
            c = min(1.0, 36 * np.pi * V ** 2 / S ** 3)
            out.update(compactness_3d=round(c, 3), sphericity=round(c ** (1 / 3), 3))
    except Exception:
        pass
    P = np.argwhere(mm) * vox
    ev = np.sort(np.linalg.eigvalsh(np.cov((P - P.mean(0)).T)))[::-1] if len(P) > 3 else None
    if ev is not None:
        a = np.sqrt(np.maximum(ev, 1e-6))
        sd_min = max(a[2], vox.min() / np.sqrt(12))                     # a 1-voxel-thin object: voxel size floor
        out.update(elongation_3d=round(float(a[0] / sd_min), 2), flatness_3d=round(float(a[1] / sd_min), 2))
    return out


_OFF = np.array([(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) if (a, b, c) > (0, 0, 0)])


def path_features(mm, vox):
    """Family E: two-sweep geodesic diameter through the candidate's voxels (largest 26-connected piece), and
    tortuosity = geodesic length / straight end-to-end distance (1 = straight; curved veins and wrapping cSS > 1)."""
    vox = np.asarray(vox, float)
    lab, n = ndi.label(mm, structure=np.ones((3, 3, 3), bool))       # merged objects: pieces need not touch;
    if n > 1: mm = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)   # the path runs in the largest piece
    idx = np.argwhere(mm)
    if len(idx) < 3: return dict(path_length_mm=np.nan, tortuosity=np.nan)
    lin = np.ravel_multi_index(idx.T, mm.shape); order = np.argsort(lin); sl = lin[order]
    r_, c_, w_ = [], [], []
    for o in _OFF:
        q = idx + o; ok = np.all((q >= 0) & (q < np.array(mm.shape)), axis=1)
        ql = np.full(len(idx), -1, np.int64); ql[ok] = np.ravel_multi_index(q[ok].T, mm.shape)
        pos = np.clip(np.searchsorted(sl, ql), 0, len(sl) - 1); hit = ok & (sl[pos] == ql)
        r_.append(np.nonzero(hit)[0]); c_.append(order[pos[hit]]); w_.append(np.full(hit.sum(), float(np.linalg.norm(o * vox))))
    r, c, w = np.concatenate(r_), np.concatenate(c_), np.concatenate(w_)
    G = coo_matrix((np.r_[w, w], (np.r_[r, c], np.r_[c, r])), shape=(len(idx),) * 2).tocsr()
    d0 = dijkstra(G, indices=0); e1 = int(np.argmax(np.where(np.isfinite(d0), d0, -1)))
    d1 = dijkstra(G, indices=e1); e2 = int(np.argmax(np.where(np.isfinite(d1), d1, -1)))
    L = float(d1[e2]); eu = float(np.linalg.norm((idx[e1] - idx[e2]) * vox))
    return dict(path_length_mm=round(L, 1), tortuosity=round(L / eu, 2) if eu > 1.0 else np.nan)


def darkness_features(zc):
    """Family A: lower tail and dark fractions of the candidate's local darkness z"""
    zc = np.asarray(zc, float)
    if not len(zc): return dict(darkness_z_p10=np.nan, dark_frac_z2=np.nan, dark_frac_z3=np.nan, dark_frac_z4=np.nan)
    return dict(darkness_z_p10=round(float(-np.percentile(zc, 10)), 2), dark_frac_z2=round(float((zc < -2).mean()), 2),
                dark_frac_z3=round(float((zc < -3).mean()), 2), dark_frac_z4=round(float((zc < -4).mean()), 2))


def sulcal_depth_map(tissue, vox, close_mm=5.0):
    """distance (mm) below the 'hull' of the brain tissue (closing by close_mm): ~0 on gyral crowns, grows into sulci"""
    d_out = ndi.distance_transform_edt(~tissue, sampling=vox)
    hull = ndi.distance_transform_edt(d_out <= close_mm, sampling=vox) > close_mm
    hull |= tissue
    return ndi.distance_transform_edt(hull, sampling=vox).astype(np.float32)


def qc_flags(row):
    """';'-joined reasons why some of this candidate's features are less trustworthy"""
    f = []
    if row.get("ref_valid_frac", 1) is not None and np.isfinite(row.get("ref_valid_frac", 1)) and row.get("ref_valid_frac", 1) < 0.5:
        f.append("cortex_reference_sparse")
    mv = row.get("mirror_valid_frac", np.nan)
    if np.isfinite(mv) and mv < 0.5: f.append("mirror_outside_brain")
    if row.get("missing_phase_features"): f.append("no_phase")
    if row.get("missing_flair_features"): f.append("no_flair")
    if row.get("missing_t1_anatomy"): f.append("swi_only_anatomy")
    if row.get("mixed_object"): f.append("mixed_merged_object")
    return ";".join(f)
