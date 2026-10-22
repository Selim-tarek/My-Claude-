"""Per-case quality control for detect_css.py (v4.14) -> review/ID_qc.json.

Called by detect_css.py with the arrays it already has. Records the acquisition, orientation, coverage,
anatomy source and how well the label (pial) boundary agrees with the cortex/CSF edge seen on the SWI,
and sets a case status:
  ok                     nothing below triggered
  low_confidence_review  anatomy-based features are less trustworthy (SWI-only anatomy, poor edge
                         agreement, high registration cost, thick slices)
  insufficient_quality   the detector's assumptions do not hold (array axis 2 is not the slice /
                         superior-inferior axis, or labels cover too little of the brain)
v4.14: REPORTED ONLY. The status and the thresholds below are provisional (not validated on expert-
labelled cases) and do not change any candidate, exclusion or score."""
import os, json, glob, numpy as np, nibabel as nib
from scipy import ndimage as ndi

EDGE_LOW_CONF_MM = 1.5     # median |label pial boundary - SWI cortex/CSF edge| above this -> low confidence
REG_COST_LOW_CONF = 0.8    # bbregister mincost above this -> low confidence (FreeSurfer typical < ~0.6)
LABEL_SLICE_MIN = 0.9      # fraction of brain-containing slices that also contain labels
THICK_SLICE_MM = 2.0       # slices thicker than this -> low confidence (sub-voxel pial thresholds)
EDGE_WIN_MM = 3.0          # search for the SWI edge within +-3 mm of the label boundary, along the normal


def read_acq(base, subj):
    """acquisition parameters from the dcm2niix JSON of the (non-phase) SWI series, if present"""
    for jp in sorted(glob.glob(f"{base}/raw/{subj}/SWI/*.json")):
        try:
            js = json.load(open(jp))
        except Exception:
            continue
        it = " ".join(js.get("ImageType", [])).upper()
        if "PHASE" in it or " P " in f" {it} ":
            continue
        te = js.get("EchoTime")
        return dict(sequence="SWI" if "SWI" in it else ("T2*-GRE" if it else "unknown"),
                    image_type=it, manufacturer=js.get("Manufacturer", ""),
                    field_T=js.get("MagneticFieldStrength"),
                    TE_ms=round(te * 1000, 1) if isinstance(te, (int, float)) else None,
                    slice_thickness_mm=js.get("SliceThickness"), acq_dim=js.get("MRAcquisitionType"))
    return {}


def edge_agreement(If, brain, pial_sd, normal, vox, sel, rng, normals=False):
    """Signed offset (mm) of the strongest dark->bright step (cortex -> CSF) along the outward normal,
    searched within +-EDGE_WIN_MM of label-boundary voxels. 0 = the labels' pial boundary sits on the
    SWI edge. Returns offsets of up to 4000 sampled boundary voxels in `sel` (and, with normals=True, the
    outward unit normals of those voxels, for the shift fit)."""
    pts = np.argwhere(sel)
    if len(pts) < 50:
        return (np.array([]), np.zeros((0, 3))) if normals else np.array([])
    if len(pts) > 4000:
        pts = pts[rng.choice(len(pts), 4000, replace=False)]
    n = np.stack([normal[a][tuple(pts.T)] for a in range(3)], 1)
    t = np.arange(-EDGE_WIN_MM, EDGE_WIN_MM + 1e-6, 0.25, dtype=np.float32)
    q = pts[:, None, :] + t[None, :, None] * n[:, None, :] / np.array(vox)
    co = q.reshape(-1, 3).T
    prof = ndi.map_coordinates(If, co, order=1, mode="nearest").reshape(q.shape[:2])
    inb = ndi.map_coordinates(brain.astype(np.float32), co, order=0, mode="constant").reshape(q.shape[:2]) > 0
    d = np.diff(prof, axis=1)
    d[~(inb[:, 1:] & inb[:, :-1])] = -np.inf          # steps that leave the brain mask do not count
    j = np.argmax(d, axis=1); ok = np.isfinite(d[np.arange(len(d)), j]) & (d[np.arange(len(d)), j] > 0)
    off = (0.5 * (t[j] + t[j + 1]))[ok]
    return (off, n[ok]) if normals else off


def _offsets(If, brain, normal, vox, pts, centre):
    """offset (mm, along the outward normal) of the strongest dark->bright step within centre +- EDGE_WIN_MM"""
    n = np.stack([normal[a][tuple(pts.T)] for a in range(3)], 1)
    t = np.arange(-EDGE_WIN_MM, EDGE_WIN_MM + 1e-6, 0.25, dtype=np.float32)
    tt = centre[:, None] + t[None, :]
    q = pts[:, None, :] + tt[:, :, None] * n[:, None, :] / np.array(vox)
    co = q.reshape(-1, 3).T
    prof = ndi.map_coordinates(If, co, order=1, mode="nearest").reshape(q.shape[:2])
    inb = ndi.map_coordinates(brain.astype(np.float32), co, order=0, mode="constant").reshape(q.shape[:2]) > 0
    d = np.diff(prof, axis=1)
    d[~(inb[:, 1:] & inb[:, :-1])] = -np.inf
    j = np.argmax(d, axis=1); r = np.arange(len(d))
    ok = np.isfinite(d[r, j]) & (d[r, j] > 0)
    off = np.where(ok, 0.5 * (tt[r, j] + tt[r, j + 1]), np.nan)
    return off, n


def shift_fit(If, brain, normal, vox, sel, rng, iters=6):
    """Rigid label->SWI shift s (mm): a label map displaced by s moves its boundary by s, so the SWI edge appears
    at -s.n along each outward normal n. Iterative: each search window is re-centred on the current estimate
    (a fixed +-3 mm window cannot see a large shift), offsets > 1.5 mm from the fit are dropped as outliers
    (veins, the white-matter edge). Returns (s, bias, scatter, share of samples with an edge): bias = median
    offset left after the shift (labels too small > 0 / too large < 0), scatter = median absolute residual."""
    pts = np.argwhere(sel)
    if len(pts) < 50:
        return None, None, None, None
    if len(pts) > 4000:
        pts = pts[rng.choice(len(pts), 4000, replace=False)]
    s = np.zeros(3)
    for k in range(iters):
        n = np.stack([normal[a][tuple(pts.T)] for a in range(3)], 1)
        off, n = _offsets(If, brain, normal, vox, pts, -(n @ s))
        ok = np.isfinite(off)
        if ok.sum() < 50:
            return None, None, None, float(ok.mean())
        res = off + n @ s
        if k:                                          # robust: drop outliers around the current fit
            ok &= np.abs(res - np.median(res[ok])) < 1.5
        s_new, *_ = np.linalg.lstsq(n[ok], -off[ok], rcond=None)
        done = np.linalg.norm(s_new - s) < 0.05
        s = s_new
        if done:
            break
    res = (off + n @ s)[np.isfinite(off)]
    return s, float(np.median(res)), float(np.median(np.abs(res - np.median(res)))), float(np.isfinite(off).mean())


def case_qc(base, subj, swi, I, If, seg, brain, tissue, cortex, pial_sd, normal, rim, hemiL, vox,
            t1_anat, full_bottom):
    rng = np.random.default_rng(0)
    qc = dict(subject=subj, qc_version="v4.19.3")
    qc["acquisition"] = read_acq(base, subj)
    qc["voxel_mm"] = [round(v, 3) for v in vox]
    qc["orientation"] = "".join(nib.aff2axcodes(swi.affine))
    si_axis = int(np.argmax(np.abs(swi.affine[2, :3])))
    qc["si_axis"] = si_axis
    qc["thickest_axis"] = int(np.argmax(vox))
    # detect_css.py runs its 2-D ridge filter, slice counts and merge gap in array axis 2 and treats
    # axis 2 as superior-inferior ("long structure"); a sagittally/coronally stored SWI breaks that
    qc["orientation_ok"] = bool(si_axis == 2 and (vox[2] >= max(vox[0], vox[1]) - 1e-3))
    qc["coverage_includes_skull_base"] = bool(full_bottom)
    if os.path.exists(f"{base}/work/{subj}_a2009s_swispace.nii.gz") and t1_anat:
        qc["anatomy"] = "T1 recon-all"
    else:
        qc["anatomy"] = "T1" if t1_anat else "SWI-only"
    reg_p = f"{base}/work/{subj}_reg.json"
    qc["registration"] = json.load(open(reg_p)) if os.path.exists(reg_p) else None
    # label coverage: brain-containing slices (axis 2) that also contain labels; labelled share of brain
    bs = brain.any(axis=(0, 1)); ls = (seg > 0).any(axis=(0, 1))
    qc["label_slice_coverage"] = round(float((bs & ls).sum() / max(int(bs.sum()), 1)), 3)
    qc["label_voxel_frac"] = round(float(((seg > 0) & brain).sum() / max(int(brain.sum()), 1)), 3)
    # pial confidence: label boundary vs SWI cortex/CSF edge, overall and by hemisphere x height third
    # v4.18.3: the voxel layers on both sides of the label boundary sit at |pial_sd| >= the smallest voxel size,
    # so 0.75 x the largest voxel selected NOTHING on isotropic scans (LEMON 0.8 mm: edge_n 0, "not measurable")
    bnd = cortex & (np.abs(pial_sd) <= max(0.75 * max(vox), 1.01 * min(vox))) & ~rim
    off = edge_agreement(If, brain, pial_sd, normal, vox, bnd, rng)
    s_, bias_, spread_, found_ = shift_fit(If, brain, normal, vox, bnd, rng)
    # v4.17.1: separate a systematic misregistration (rigid shift, global bias) from scatter. On P006/P011 the
    # median |offset| (1.6 / 2.1 mm) mixed both; the status now uses the systematic parts only.
    qc["label_shift_mm"] = round(float(np.linalg.norm(s_)), 2) if s_ is not None else None
    qc["label_shift_vec_mm"] = [round(float(v), 2) for v in s_] if s_ is not None else None
    qc["edge_bias_mm"] = round(bias_, 2) if bias_ is not None else None
    qc["edge_scatter_mm"] = round(spread_, 2) if spread_ is not None else None
    qc["edge_found_frac"] = round(found_, 2) if found_ is not None else None
    qc["edge_offset_median_mm"] = round(float(np.median(off)), 2) if len(off) else None
    qc["edge_agreement_mm"] = round(float(np.median(np.abs(off))), 2) if len(off) else None
    qc["edge_n"] = int(len(off))
    zw = (swi.affine @ np.c_[np.argwhere(bnd), np.ones(int(bnd.sum()))].T)[2] if bnd.any() else np.array([])
    reg = {}
    if len(zw):
        lo, hi = np.percentile(zw, [0, 100]); cuts = np.linspace(lo, hi, 4)
        zmap = np.full(I.shape, -1, np.int8); zmap[bnd] = np.clip(np.searchsorted(cuts, zw, side="right") - 1, 0, 2)
        for h, hm in (("L", hemiL), ("R", ~hemiL)):
            for k, nm in enumerate(("inferior", "middle", "superior")):
                o = edge_agreement(If, brain, pial_sd, normal, vox, bnd & hm & (zmap == k), rng)
                reg[f"{h}_{nm}"] = round(float(np.median(np.abs(o))), 2) if len(o) >= 50 else None
    qc["edge_agreement_by_region_mm"] = reg
    # motion / ghosting proxy (informational): background noise of the un-stripped SWI vs brain signal
    raw_p = f"{base}/work/{subj}_raw.nii.gz"
    qc["background_ratio"] = None
    if os.path.exists(raw_p):
        raw = nib.load(raw_p)
        if raw.shape[:3] == I.shape:
            R = np.asarray(raw.dataobj, dtype=np.float32).reshape(I.shape)
            far = ndi.distance_transform_edt(~brain, sampling=vox) > 15.0
            if far.sum() > 1000 and brain.any():
                qc["background_ratio"] = round(float(np.std(R[far]) / (np.median(R[brain]) + 1e-6)), 4)
    # status (provisional thresholds, reported only)
    why_bad, why_low = [], []
    if not qc["orientation_ok"]:
        why_bad.append(f"array axis 2 is not the slice / S-I axis (orientation {qc['orientation']}, voxel {qc['voxel_mm']})")
    if qc["label_slice_coverage"] < LABEL_SLICE_MIN:
        why_bad.append(f"labels cover {qc['label_slice_coverage']:.0%} of brain slices")
    if qc["anatomy"] == "SWI-only":
        why_low.append("SWI-only anatomy (no T1)")
    if qc["edge_agreement_mm"] is None:
        why_low.append("pial edge agreement not measurable")
    else:
        if qc["label_shift_mm"] is not None and qc["label_shift_mm"] > EDGE_LOW_CONF_MM:
            why_low.append(f"labels shifted {qc['label_shift_mm']} mm against the SWI (rigid fit to the cortex edge)")
        if qc["edge_bias_mm"] is not None and abs(qc["edge_bias_mm"]) > EDGE_LOW_CONF_MM:
            why_low.append(f"label pial boundary {qc['edge_bias_mm']:+.2f} mm from the SWI cortex edge on average "
                           f"({'labels too small' if qc['edge_bias_mm'] > 0 else 'labels too large'})")
    cost = (qc["registration"] or {}).get("cost")
    # v4.19.3: the bbregister cost only lowers the status when the direct edge fit could not be measured. LEMON GRE
    # (P201/P202/P204): cost 0.95-0.96 every time while the labels sat 0.11-0.22 mm from the SWI cortex edge -
    # the T2*-weighted GRE contrast gives high costs even when aligned. The cost is still recorded.
    qc["notes"] = []
    if isinstance(cost, (int, float)) and cost > REG_COST_LOW_CONF:
        if qc["label_shift_mm"] is None:
            why_low.append(f"registration cost {cost}")
        else:
            qc["notes"].append(f"registration cost {cost} above {REG_COST_LOW_CONF} - not used, the edge fit was measured "
                               f"(label shift {qc['label_shift_mm']} mm)")
    if max(vox) > THICK_SLICE_MM:
        why_low.append(f"thick slices ({max(vox):.1f} mm)")
    qc["status"] = "insufficient_quality" if why_bad else ("low_confidence_review" if why_low else "ok")
    qc["reasons"] = why_bad + why_low
    qc["thresholds"] = dict(edge_low_conf_mm=EDGE_LOW_CONF_MM, reg_cost_low_conf=REG_COST_LOW_CONF,
                            label_slice_min=LABEL_SLICE_MIN, thick_slice_mm=THICK_SLICE_MM,
                            note="provisional, not validated; v4.14 reports only")
    os.makedirs(f"{base}/review", exist_ok=True)
    json.dump(qc, open(f"{base}/review/{subj}_qc.json", "w"), indent=2)
    return qc
