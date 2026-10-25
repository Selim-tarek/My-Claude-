"""Review-priority tiers (v5.1.3): a transparent READING ORDER from the recorded evidence. Not a diagnosis, not an
exclusion: every candidate stays in the list and in the audit trail. Pure function; unit-tested in tests/check_v513.py.

  1 high_priority_cSS_candidate   no strong mimic evidence
  2 css_vs_mimic_uncertain        mixed merged object; possible calcification (phase); FLAIR context (acute cSAH /
                                  thrombosis)
  3 likely_vein                   CONVERGENT vascular evidence (never one feature): white-matter continuation >= 4 mm
                                  with dark continuation >= 8 mm or a tubular shape (transcortical); reaches the sinus
                                  zone with >= 10 mm continuation; tubular and continuing at both ends; OFF the cortex
                                  (pial distance > 0.3 mm) AND in the CSF (>= 0.5 mm from tissue) AND vessel-shaped
                                  (tube ratio > 0.35, elongation >= 10 or continuation >= 10 mm). cSS in the
                                  subarachnoid space can sit off the cortex too - one reason this only reorders.
  4 likely_other_mimic_or_artifact  air-bone artifact zone; microbleed-like (single compact parenchymal focus or cluster)
  5 case_qc_insufficient          the case failed QC (insufficient_quality): fix the case before reading
A low_confidence_review case keeps its per-candidate tiers (a case-wide constant would not order anything); the
case status is in qc_status and the reviewer shows it.
All cut-offs are PLACEHOLDERS from the v4 rules and phantoms, not validated; compare the order with expert labels
(evaluate.py) before relying on it."""
import numpy as np

TIER_NAMES = {1: "high_priority_cSS_candidate", 2: "css_vs_mimic_uncertain", 3: "likely_vein",
              4: "likely_other_mimic_or_artifact", 5: "case_qc_insufficient"}
WM_CONT_MM, WM_RUN_MM, SINUS_RUN_MM, TUBE_HI, TUBE_MID, CSF_DEPTH_MM, ELONG_HI = 4.0, 8.0, 10.0, 0.45, 0.35, 0.5, 10.0


def _f(r, k, d=np.nan):
    v = r.get(k, d)
    try: return float(v)
    except (TypeError, ValueError): return d


def review_priority(r, qc_status):
    """candidate row (dict) -> (tier 1-5, name, reason)"""
    if qc_status == "insufficient_quality":
        return 5, TIER_NAMES[5], "case QC insufficient_quality"
    if _f(r, "mixed_object", 0) >= 1:
        return 2, TIER_NAMES[2], "mixed merged object (vein-like and cSS-like pieces)"
    if str(r.get("mimic_class", "none")) in ("mb_cluster",) or _f(r, "cmb_like", 0) >= 1:
        return 4, TIER_NAMES[4], "microbleed-like (compact parenchymal focus / cluster)"
    if _f(r, "artifact_zone", 0) >= 1:
        return 4, TIER_NAMES[4], "air-bone artifact zone"
    tube, cont = _f(r, "tube_ratio"), _f(r, "local_cont_max_mm", 0)
    vein = []
    if _f(r, "vessel_wm_mm", 0) >= WM_CONT_MM and (cont >= WM_RUN_MM or (np.isfinite(tube) and tube > TUBE_HI)):
        vein.append("continues into white matter")
    if "midline" in (str(r.get("end_type_max", "")), str(r.get("end_type_min", ""))) and cont >= SINUS_RUN_MM:
        vein.append("runs on to the sinus zone")
    if np.isfinite(tube) and tube > TUBE_HI and _f(r, "terminal_ends", 2) == 0: vein.append("tubular, continues at both ends")
    shaped = (np.isfinite(tube) and tube > TUBE_MID) or _f(r, "elongation_3d", 0) >= ELONG_HI or cont >= SINUS_RUN_MM
    if _f(r, "pial_dist_mm", 0) > 0.3 and _f(r, "csf_depth_mm", 0) >= CSF_DEPTH_MM and shaped:
        vein.append("off the cortex, in the CSF, vessel-shaped")
    if str(r.get("mimic_class", "none")) == "calcification":
        return 2, TIER_NAMES[2], "possible calcification (diamagnetic phase)" + (f"; also {vein[0]}" if vein else "")
    if str(r.get("mimic_class", "none")) in ("acute", "thrombosis", "mixed"):
        return 2, TIER_NAMES[2], f"FLAIR context ({r.get('mimic_class')})"
    if vein:
        return 3, TIER_NAMES[3], "; ".join(vein)
    return 1, TIER_NAMES[1], "no strong mimic evidence"
