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
