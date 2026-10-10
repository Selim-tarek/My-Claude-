"""Shared helpers for the cSS pipeline scripts (imported from the same folder)."""
import os, numpy as np, pandas as pd

FP = ["vox_i", "vox_j", "vox_k", "volume_mm3"]      # candidate fingerprint stored with each call


def base():
    return os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))


def load_lut():
    lut = {}
    p = os.path.join(os.environ.get("FREESURFER_HOME", ""), "FreeSurferColorLUT.txt")
    if os.path.exists(p):
        for line in open(p):
            t = line.split()
            if len(t) > 2 and t[0].isdigit():
                lut[int(t[0])] = t[1]
    return lut


def match_calls(cand, calls, tol_vox=2.0, tol_vol=0.15):
    """Re-attach reader calls to the CURRENT candidate list by fingerprint (centroid within
    tol_vox voxels, volume within tol_vol), so re-running the detector can never move a call
    onto a different lesion. Returns ({cand_id: call}, n_dropped), or (None, n) for a legacy
    calls file that has no fingerprint columns."""
    if not all(c in calls.columns for c in FP):
        return None, len(calls)
    out, dropped = {}, 0
    C = cand[["cand_id"] + FP].to_numpy(float)
    for _, r in calls.iterrows():
        d = np.sqrt(((C[:, 1:4] - r[FP[:3]].to_numpy(float)) ** 2).sum(1))
        dv = np.abs(C[:, 4] - float(r.volume_mm3)) / max(float(r.volume_mm3), 1.0)
        ok = np.where((d <= tol_vox) & (dv <= tol_vol))[0]
        if len(ok):
            out[int(C[ok[d[ok].argmin()], 0])] = r.call
        else:
            dropped += 1
    return out, dropped


def load_calls(b, subj, cand):
    """Reader calls for subj mapped onto the current candidates. Returns (dict, warning)."""
    p = f"{b}/review/{subj}_calls.csv"
    if not os.path.exists(p) or not len(cand):
        return {}, ""
    calls = pd.read_csv(p)
    if not len(calls):
        return {}, ""
    m, dropped = match_calls(cand, calls)
    if m is None:
        # v3 file (cand_id only): trust it only if the detector has not been re-run since.
        # candidates.nii.gz is written by detect_css.py only (mark_css.py rewrites the csv).
        if os.path.getmtime(f"{b}/review/{subj}_candidates.nii.gz") > os.path.getmtime(p):
            return {}, (f"{subj}: old-format calls IGNORED - candidates were re-detected after "
                        f"they were made, so cand_ids may point to different lesions")
        valid = set(cand.cand_id.astype(int))
        return {int(c): v for c, v in zip(calls.cand_id, calls.call) if int(c) in valid}, ""
    return m, (f"{subj}: {dropped} earlier call(s) match no current candidate (detector re-run) "
               f"- dropped" if dropped else "")


def save_calls(b, subj, cand, calls):
    c = cand.set_index("cand_id")
    rows = [dict(cand_id=k, call=v, **{f: c.loc[k, f] for f in FP}) for k, v in sorted(calls.items())]
    pd.DataFrame(rows, columns=["cand_id", "call"] + FP).to_csv(f"{b}/review/{subj}_calls.csv", index=False)


# ---- v4.15 expert label vocabulary (docs/AUDIT_v4.13_and_FP_plan.md section 4)
# sheet code -> (fine label, coarse call). The coarse call is what review_css.py / feature_report.py /
# mark_css.py already understand, so old single letters (C V N H U) keep their meaning.
LABELS = {
    "C":   ("cSS", "cSS"),
    "V":   ("vein_unspecified", "Vein"),
    "VS":  ("vein_surface", "Vein"),            # cortical vein running ALONG the pial surface (hard negative)
    "VC":  ("vein_sulcal", "Vein"),             # vein in the middle of the sulcal CSF
    "VT":  ("vein_transcortical", "Vein"),      # medullary / transcortical, perpendicular to the cortex
    "VX":  ("vein_sinus_adjacent", "Vein"),     # vein joining / next to a dural sinus
    "TV":  ("thrombosed_vein", "Vein"),
    "A":   ("artifact_airbone", "Artifact"),    # skull base / orbitofrontal / temporal pole susceptibility
    "AM":  ("artifact_motion", "Artifact"),     # motion, ghosting, Gibbs, registration
    "N":   ("normal_other", "Normal"),
    "NI":  ("normal_iron_cortex", "Normal"),    # e.g. iron-rich motor cortex
    "MB":  ("microbleed", "Normal"),
    "CA":  ("calcification", "Normal"),
    "SAH": ("acute_cSAH", "Normal"),
    "LN":  ("laminar_necrosis", "Normal"),
    "HI":  ("haemorrhagic_infarct", "Normal"),
    "IS":  ("infratentorial_SS", "Normal"),
    "H":   ("ICH_related", "Near ICH"),
    "U":   ("uncertain", "Unsure"),
}
LABEL_HELP = [("C", "cortical superficial siderosis"),
              ("V / VS / VC / VT / VX", "vein: any / along surface / mid-sulcus / transcortical / at sinus"),
              ("TV", "thrombosed vein"),
              ("A / AM", "artifact: air-bone susceptibility / motion, ghosting"),
              ("N / NI", "normal or other / normal iron-rich cortex"),
              ("MB / CA", "microbleed / calcification"),
              ("SAH", "acute convexity SAH (FLAIR-bright sulcus)"),
              ("LN / HI", "laminar necrosis / haemorrhagic infarct"),
              ("IS", "infratentorial (classical) superficial siderosis"),
              ("H", "siderosis contiguous with a lobar ICH"),
              ("U", "uncertain")]


def parse_answer(tok):
    """'12VS' or '12VS:4' -> (12, code, confidence or None); raises ValueError"""
    import re
    m = re.fullmatch(r"(\d+)([A-Z]+)(?::([1-5]))?", tok.strip().upper())
    if not m or m.group(2) not in LABELS:
        raise ValueError(f"cannot read '{tok}' (e.g. 12C, 12VS, 12VS:4; codes {', '.join(LABELS)})")
    return int(m.group(1)), m.group(2), (int(m.group(3)) if m.group(3) else None)


def cand_uid(subj, r):
    """stable candidate identity across detector re-runs (same fingerprint as match_calls)"""
    return f"{subj}:{int(r['vox_i'])}:{int(r['vox_j'])}:{int(r['vox_k'])}:{float(r['volume_mm3']):g}"


def append_rows(path, rows, key):
    """append rows to a CSV table, replacing earlier rows with the same key columns"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    new = pd.DataFrame(rows)
    if os.path.exists(path):
        old = pd.read_csv(path)
        if len(new):
            k_old = old[key].astype(str).agg("|".join, axis=1)
            k_new = set(new[key].astype(str).agg("|".join, axis=1))
            old = old[~k_old.isin(k_new)]
        new = pd.concat([old, new], ignore_index=True)
    new.to_csv(path, index=False)
    return new
