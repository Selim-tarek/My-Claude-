"""Per-piece evidence for merged candidates (v5.1.1). Pure functions, imported by detect_css.py; unit-tested in
tests/check_v511.py.

detect_css.py merges dark pieces closer than MERGE_GAP_MM into one candidate (a line broken by noise, the two banks of
a tram-track). That can also join a cSS segment and a neighbouring vein; the group's median features then describe
neither, and a vein rule (E5 / E6 / E7) would remove the cSS segment with it. Each pre-merge piece is therefore
classified with the detector's own vein_like definition:
  vein-like piece   median pial distance > 0.3 mm (off the cortex) and tubular (median tube ratio > 0.4, or unknown)
  cSS-like piece    median pial distance <= 0.3 mm and not tubular (tube ratio <= 0.4, or unknown)
mixed_object = >= 2 pieces and each kind holds >= min_frac of the object's volume.
Limits: a vein TOUCHING a cSS segment voxel to voxel is one piece already (this cannot split it); merges undone by
the 40 mm chain cap leave single pieces (nothing to compare). Thresholds are the v4 vein_like ones, not validated."""
import numpy as np


def piece_evidence(mm, pieces, pial_sd, tube, min_frac=0.15):
    """mm: candidate mask (crop); pieces: pre-merge labels (same crop, or None = no merging); pial_sd, tube: maps (same
    crop). -> dict of the per-piece columns"""
    pcs = {}
    if pieces is not None:
        pl = pieces[mm]
        for pid in np.unique(pl[pl > 0]):
            pm = mm & (pieces == pid)
            tv = tube[pm]; tp = float(np.nanmedian(tv)) if np.isfinite(tv).any() else np.nan
            pp = float(np.median(pial_sd[pm]))
            pcs[int(pid)] = (int(pm.sum()), pp, tp, bool(pp > 0.3 and (np.isnan(tp) or tp > 0.4)),
                             bool(pp <= 0.3 and (np.isnan(tp) or tp <= 0.4)))
    tot = sum(v[0] for v in pcs.values()) or 1
    vf = sum(v[0] for v in pcs.values() if v[3]) / tot
    cf = sum(v[0] for v in pcs.values() if v[4]) / tot
    return dict(piece_vol_vox=";".join(str(v[0]) for v in pcs.values()),
                piece_pial_mm=";".join(f"{v[1]:.2f}" for v in pcs.values()),
                piece_tube=";".join("nan" if np.isnan(v[2]) else f"{v[2]:.2f}" for v in pcs.values()),
                piece_vein_frac=round(vf, 2) if pcs else np.nan,
                mixed_object=int(len(pcs) > 1 and vf >= min_frac and cf >= min_frac))


def rescue_mixed(rule_version, mixed_object, enabled):
    """True = a vein rule (E5 / E6 / E7) must not remove this mixed object as a whole (CSS_MIXED_REVIEW=1)"""
    return bool(enabled and mixed_object and rule_version.split("@")[0] in ("E5", "E6", "E7"))
