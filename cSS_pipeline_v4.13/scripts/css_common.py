"""Shared helpers for the cSS pipeline scripts (imported from the same folder)."""
import os, sys, numpy as np, pandas as pd

LABEL_SCHEMA = "css-labels/1.0"   # v5.1.1: expert label table schema (fine code + label + coarse class + confidence)
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


# ---- v4.17 provenance helpers
def stale_against(path, sources):
    """v5.1: the existing sources that are NEWER than path (empty list = path is fresh; path missing = all of them)"""
    if not os.path.exists(path): return [s for s in sources if os.path.exists(s)]
    t = os.path.getmtime(path)
    return [s for s in sources if os.path.exists(s) and os.path.getmtime(s) > t]


def use_derived(path, sources, what):
    """v5.1 freshness gate for a derived file: True = use it. Stale -> ignored with a message, or (CSS_QC_STRICT=1)
    stop. Missing -> False, silently (callers decide whether it is required)."""
    if not os.path.exists(path): return False
    newer = stale_against(path, sources)
    if not newer: return True
    msg = f"{what} ({os.path.basename(path)}) is older than {', '.join(os.path.basename(s) for s in newer)}"
    if os.environ.get("CSS_QC_STRICT") == "1": sys.exit(f"ERROR: {msg} - CSS_QC_STRICT=1 stops (re-run the step that makes it)")
    print(f"warning: {msg} - ignored"); return False


def fresh(path, ref):
    """True if path exists and is not older than ref (a derived file made after its source was imported)"""
    return os.path.exists(path) and (not os.path.exists(ref) or os.path.getmtime(path) >= os.path.getmtime(ref))


def sha256(path, block=1 << 20):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


MANIFEST_SCHEMA = "css-manifest/1.0"


def build_id():
    """v5.1: the package build. git archive fills scripts/BUILD (export-subst); a source checkout asks git"""
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BUILD")
    if os.path.exists(p):
        t = open(p).read().strip()
        if t and "$Format" not in t: return t
    try:
        import subprocess
        r = subprocess.run(["git", "-C", os.path.dirname(os.path.abspath(__file__)), "rev-parse", "HEAD"],
                           capture_output=True, text=True, timeout=5)
        if r.returncode == 0: return r.stdout.strip() + " (source checkout)"
    except Exception:
        pass
    return "unknown"


def update_manifest(b, subj, section, data):
    """merge {section: data} into review/ID_manifest.json (one machine-readable provenance record per case)"""
    import json
    p = f"{b}/review/{subj}_manifest.json"
    m = json.load(open(p)) if os.path.exists(p) else {"subject": subj}
    m["schema_version"] = MANIFEST_SCHEMA            # v5.1: machine-readable schema id (docs/CHANGELOG_v5.1.md)
    m[section] = data
    os.makedirs(os.path.dirname(p), exist_ok=True)
    json.dump(m, open(p, "w"), indent=2, default=str)


def environment_record():
    """python / package / FreeSurfer versions for the manifest"""
    import platform, importlib
    rec = {"python": platform.python_version(), "platform": platform.platform(), "build": build_id()}
    for mod in ("numpy", "scipy", "nibabel", "pandas", "skimage", "matplotlib"):
        try:
            rec[mod] = importlib.import_module(mod).__version__
        except Exception:
            rec[mod] = None
    fs = os.environ.get("FREESURFER_HOME", "")
    stamp = os.path.join(fs, "build-stamp.txt")
    rec["freesurfer"] = open(stamp).read().strip() if os.path.exists(stamp) else (fs or None)
    return rec


def label_class(label):
    """fine expert label -> display / agreement class (shared by surface_overlay.py and consensus.py)"""
    if label == "cSS": return "cSS"
    if label == "uncertain": return "uncertain"
    if label == "vein_surface": return "surface_vein"
    if label.startswith("vein_") or label == "thrombosed_vein": return "other_vein"
    if label.startswith("artifact"): return "artifact"
    return "other_normal"


def cohen_kappa(a, b, categories=None, weights=None):
    """Cohen's kappa for two raters; weights=None (unweighted), 'linear' or 'quadratic' (ordered categories)"""
    a, b = list(a), list(b)
    cats = list(categories) if categories is not None else sorted(set(a) | set(b))
    if len(a) == 0 or len(cats) < 2:
        return float("nan")
    k = {c: i for i, c in enumerate(cats)}; n = len(cats)
    O = np.zeros((n, n))
    for x, y in zip(a, b):
        O[k[x], k[y]] += 1
    O /= O.sum()
    E = np.outer(O.sum(1), O.sum(0))
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    W = np.where(i == j, 0.0, 1.0) if weights is None else (np.abs(i - j) / (n - 1)) ** (2 if weights == "quadratic" else 1)
    den = (W * E).sum()
    return float(1 - (W * O).sum() / den) if den > 0 else float("nan")


def skeleton_length(sk, vox):
    """v4.18: centre-line length (mm) of a 2-D skeleton from its steps (detect_css used pixels x vox[0], which
    undercounts diagonals): orthogonal steps vox[0] / vox[1], diagonal steps their hypotenuse; a diagonal is
    skipped where an orthogonal corner already joins the two pixels"""
    import numpy as np
    sk = np.asarray(sk, bool); a, b = float(vox[0]), float(vox[1]); d = float(np.hypot(a, b))
    h = sk[1:, :] & sk[:-1, :]; v = sk[:, 1:] & sk[:, :-1]
    d1 = sk[1:, 1:] & sk[:-1, :-1] & ~(sk[1:, :-1] | sk[:-1, 1:])
    d2 = sk[1:, :-1] & sk[:-1, 1:] & ~(sk[1:, 1:] | sk[:-1, :-1])
    return float(h.sum() * a + v.sum() * b + (d1.sum() + d2.sum()) * d)


def radius_profile(mask, vox, sample=None):
    """v4.18: in-plane half-width (mm) along the 2-D centre line, slice by slice (cf. Avadiappan 2020: radius from a
    Euclidean distance transform). Distance to the nearest outside voxel (sampling vox[0], vox[1]) minus half a voxel,
    read on the skeleton of `sample` (default: the mask itself), so a 1-voxel-wide line reads half a voxel and a
    3-voxel-wide band 1.5 voxels. Quantised at the in-plane voxel size: compare only within one protocol."""
    import numpy as np
    from scipy import ndimage as ndi
    from skimage.morphology import skeletonize
    mask = np.asarray(mask, bool); sample = mask if sample is None else np.asarray(sample, bool)
    half = 0.5 * float(np.sqrt(vox[0] * vox[1])); out = []
    for k in range(mask.shape[2]):
        s = sample[:, :, k] & mask[:, :, k]
        if not s.any(): continue
        e = ndi.distance_transform_edt(np.pad(mask[:, :, k], 1), sampling=vox[:2])[1:-1, 1:-1]
        sk = skeletonize(s) & s
        out.append(e[sk] - half)
    return np.concatenate(out) if out else np.zeros(0)


GAMMA_HZ_T = 42.577e6       # proton gyromagnetic ratio / 2 pi


def susceptibility_signal(chi_ppm, vox, affine, b0_t=3.0, te_s=0.020, pad=16):
    """v4.18 susceptibility forward model for synthetic lesions (cf. CenSynCMB 2026: dipole kernel + static dephasing).
    chi_ppm: susceptibility DIFFERENCE map (ppm, 0 = host tissue). The field perturbation (ppm of B0) is the
    convolution with the unit dipole, done in k-space: D(k) = 1/3 - (k.b0)^2 / |k|^2, with b0 the scanner z axis
    expressed in the array axes (oblique acquisitions handled through the affine). Phase = 2 pi gamma B0 TE field.
    Intravoxel dephasing for a field that varies linearly across the voxel: |sinc| along each axis of the local
    phase gradient x voxel size. Returns (phase_rad, dephasing_attenuation 0-1). An approximation for
    augmentation and stress tests - not a full GRE / SWI / QSM simulator."""
    import numpy as np
    chi = np.pad(np.asarray(chi_ppm, np.float64), pad)
    vox = np.asarray(vox, float)
    R = np.asarray(affine, float)[:3, :3] / vox[None, :]
    b0 = R.T @ np.array([0.0, 0.0, 1.0]); b0 /= np.linalg.norm(b0)          # B0 (world z) in array axes
    k = np.meshgrid(*[np.fft.fftfreq(n, d=v) for n, v in zip(chi.shape, vox)], indexing="ij")
    k2 = k[0] ** 2 + k[1] ** 2 + k[2] ** 2
    kb = k[0] * b0[0] + k[1] * b0[1] + k[2] * b0[2]
    with np.errstate(invalid="ignore", divide="ignore"):
        D = 1.0 / 3.0 - np.where(k2 > 0, kb ** 2 / k2, 0.0)
    D[0, 0, 0] = 0.0
    field = np.real(np.fft.ifftn(np.fft.fftn(chi) * D))[pad:-pad, pad:-pad, pad:-pad]       # ppm
    phase = 2 * np.pi * GAMMA_HZ_T * b0_t * te_s * field * 1e-6
    att = np.ones_like(phase)
    for ax in range(3):
        g = np.gradient(phase, vox[ax], axis=ax)                                # rad / mm
        att *= np.abs(np.sinc(g * vox[ax] / (2 * np.pi)))
    return phase.astype(np.float32), att.astype(np.float32)


# ---- v5.0 phase helpers (moved unchanged from make_swi.py so the detector's mimic layer uses the same filter)
def phase_to_radians(P):
    """phase image -> (radians float32, description). Scanner units (range wider than 2 pi + 0.5, e.g.
    -4096..4095 or 0..4095) are rescaled linearly to -pi..pi; anything else is taken to be radians already."""
    import numpy as np
    P = np.asarray(P, np.float32)
    lo, hi = float(np.nanmin(P)), float(np.nanmax(P))
    if hi - lo > 2 * np.pi + 0.5:
        return (P - lo) / (hi - lo) * 2 * np.pi - np.pi, f"scanner units {lo:g}..{hi:g} -> radians"
    return P, "radians"


def homodyne_highpass(P_rad, vox, si, sigma_mm=4.0):
    """wrap-safe homodyne high-pass, slice by slice in the plane perpendicular to array axis `si` (the S-I axis):
    angle(z * conj(lowpass z)), Gaussian low-pass sigma_mm in-plane, 0 along si (as in make_swi.py since v4.18)"""
    import numpy as np
    from scipy import ndimage as ndi
    z = np.exp(1j * np.asarray(P_rad).astype(np.float64))
    sig = [0.0, 0.0, 0.0]
    for ax in range(3):
        if ax != si: sig[ax] = sigma_mm / vox[ax]
    lp = ndi.gaussian_filter(z.real, sig) + 1j * ndi.gaussian_filter(z.imag, sig)
    return np.angle(z * np.conj(lp)).astype(np.float32)

