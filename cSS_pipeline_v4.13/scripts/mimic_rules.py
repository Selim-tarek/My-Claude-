"""Mimic evidence for detect_css.py (v5.0, Layer 1). Imported, not run on its own.

Three published cSS mimics, each measured per candidate (Charidimou AJNR 2016; mimics review, Clin Neuroradiol 2020):

  calcification (Sturge-Weber, calcifying angiopathy)   phase: DIAMAGNETIC; cSS / veins are paramagnetic
  acute convexity SAH / cortical vein thrombosis        FLAIR: bright sulcal CSF / bright adjacent cortex
  cluster of lobar microbleeds                          shape: several compact parenchymal blobs merged into one object

Columns (detect_css.py appends them to candidates AND excluded):
  mimic_score   0-1, the strongest of the three pieces of evidence below (0 = none); NOT a probability
  mimic_class   none / calcification / acute / thrombosis / mixed / mb_cluster
  mimic_reason  human-readable evidence
  phase_z       approximate t-statistic of the candidate's fitted susceptibility from the phase (> 0 paramagnetic,
                < 0 diamagnetic); NaN without a usable phase image
  phase_r       correlation of the observed filtered phase with the dipole pattern a PARAMAGNETIC object of this
                candidate's own shape would make (-1..1)
  mb_cluster_n  compact parenchymal sub-objects left after one erosion (>= 3 -> microbleed cluster)

WHY THE PHASE RULE FITS A DIPOLE PATTERN instead of reading the sign of the median phase: the phase inside a thin
sheet or tube depends on its orientation to B0. In css_common.susceptibility_signal (the pipeline's own forward
model), a paramagnetic sheet lying flat to B0 has median filtered phase -0.94 rad and a vertical one +2.09 rad; a
paramagnetic tube perpendicular to B0 -0.54 rad (tests/check_mimic.py). A sign-of-median rule would call real cSS on
horizontal cortex "calcium". Fitting the pattern a paramagnetic object of the candidate's own shape would produce
takes the orientation into account. The scanner's sign convention is fitted the same way on the vessel map (veins are
paramagnetic), so no vendor table is needed.

Limits (research software, thresholds NOT validated on labelled data):
  - phase_z treats noise as independent within ~2 voxels; real phase noise and nearby veins make it optimistic.
  - background fields that survive the 4 mm high-pass are noise here. Calcium next to a large vein can be missed.
  - FLAIR thresholds (z > 3) are the detector's existing flair_bright / flair_ctx_bright cut-offs (v4.x), not
    re-derived. Acute cSAH and cSS co-exist in CAA: FLAIR evidence is a SOFT penalty only, never an exclusion.
  - laminar cortical necrosis and haemorrhagic infarct are not handled (no DWI / T1 hyperintensity feature).
Only with CSS_USE_MIMIC=1 does any of this change the output: score x 0.5 for FLAIR evidence, exclusion E10
(calcification, mimic_score >= 0.7) and E11 (microbleed cluster, mimic_score >= 0.5). Default: recorded only."""
import os, numpy as np, nibabel as nib
from scipy import ndimage as ndi
from css_common import phase_to_radians, homodyne_highpass, fresh

HP_SIGMA_MM = 4.0          # same high-pass as make_swi.py
FIT_SHELL_VOX = 2          # phase fit region = this many voxels around the candidate (inside the brain)
WRAP_RAD = 2.5             # observed |filtered phase| above this may be wrapped: not used in the fit
BOX_MARGIN_MM = 12.0       # forward-model box around the candidate
CORR_VOX = 2.0             # assumed noise correlation length (voxels per axis) for phase_z
CONV_MIN_T = 3.0           # scanner sign convention must be fitted at least this clearly, else the phase rule is off
CALC_Z, CALC_R = -3.0, -0.2
FLAIR_Z = 3.0              # = flair_bright / flair_ctx_bright in detect_css.py
MB_DIAM_MM, MB_MIN_VOX, MB_MAX_ASPECT, MB_PARENCHYMA = 6.0, 2, 2.5, 0.5
MB_MIN_N = 3


def dipole_field(chi, vox, affine, pad=16):
    """field perturbation (ppm of B0) of susceptibility map chi: the same unit-dipole kernel as
    css_common.susceptibility_signal, D(k) = 1/3 - (k.b0)^2/|k|^2, B0 = world z through the affine; lean memory
    (broadcast frequencies, complex64) so a whole-brain fit fits on a laptop"""
    c = np.pad(np.asarray(chi, np.float32), pad)
    vox = np.asarray(vox, float)
    R = np.asarray(affine, float)[:3, :3] / vox[None, :]
    b0 = R.T @ np.array([0.0, 0.0, 1.0]); b0 /= np.linalg.norm(b0)
    f = [np.fft.fftfreq(n, d=v).astype(np.float32) for n, v in zip(c.shape, vox)]
    kx, ky, kz = f[0][:, None, None], f[1][None, :, None], f[2][None, None, :]
    k2 = kx ** 2 + ky ** 2 + kz ** 2
    kb = kx * b0[0] + ky * b0[1] + kz * b0[2]
    with np.errstate(invalid="ignore", divide="ignore"):
        D = (np.float32(1 / 3) - np.where(k2 > 0, kb ** 2 / k2, 0)).astype(np.float32)
    D[0, 0, 0] = 0
    from scipy import fft as sfft                                  # keeps complex64 (numpy.fft upcasts)
    out = np.real(sfft.ifftn(sfft.fftn(c.astype(np.complex64), workers=-1) * D, workers=-1)).astype(np.float32)
    return out[pad:-pad, pad:-pad, pad:-pad]


def highpass_linear(field, vox, si, sigma_mm=HP_SIGMA_MM):
    """the linear (small-phase) equivalent of homodyne_highpass, applied to a predicted field"""
    sig = [0.0 if ax == si else sigma_mm / vox[ax] for ax in range(3)]
    return field - ndi.gaussian_filter(field, sig)


def fit_amplitude(obs, pred, n_corr=CORR_VOX ** 3):
    """least-squares a in obs = a * pred over the given voxels -> (a, t, r); t uses n / n_corr effective samples"""
    obs = np.asarray(obs, float); pred = np.asarray(pred, float)
    ss = float((pred ** 2).sum())
    if len(obs) < 10 or ss <= 0: return np.nan, np.nan, np.nan
    a = float((obs * pred).sum() / ss)
    res = obs - a * pred
    # var(a) = sigma^2 / sum(pred^2) for independent noise; noise correlated over n_corr voxels inflates it n_corr-fold
    se = float(np.sqrt((res ** 2).sum() / max(len(obs) - 1, 1)) * np.sqrt(n_corr / ss))
    r = float(np.corrcoef(obs, pred)[0, 1]) if obs.std() > 0 and pred.std() > 0 else np.nan
    return a, (a / se if se > 0 else np.nan), r


def flair_rule(fz, fcz):
    """FLAIR evidence -> (score 0-1, class, reason). fz = sulcal-CSF FLAIR z, fcz = adjacent-cortex FLAIR z
    (detect_css.py flair_csf_z / flair_ctx_z; NaN without FLAIR)"""
    s = lambda z: float(np.clip(0.5 + (z - FLAIR_Z) / 6.0, 0, 1)) if np.isfinite(z) and z > FLAIR_Z else 0.0
    sc, sx = s(fz), s(fcz)
    if sc and sx: return max(sc, sx), "mixed", f"FLAIR bright sulcal CSF (z {fz:.1f}) and cortex (z {fcz:.1f})"
    if sc: return sc, "acute", f"FLAIR bright sulcal CSF (z {fz:.1f}): acute convexity SAH?"
    if sx: return sx, "thrombosis", f"FLAIR bright adjacent cortex (z {fcz:.1f}): cortical vein thrombosis?"
    return 0.0, "none", ""


def mb_cluster_rule(mask, tissue, vox):
    """candidate support (crop) -> (score, n, reason). Erode once with the 6-neighbour cross (a full 3x3x3 cube
    shrinks a 3.2 mm microbleed at 0.8 mm voxels to one voxel), in-plane only when slices are > 1.5x the in-plane
    voxel (at 0.5 x 0.5 x 2 mm a 3D erosion erases every microbleed under 5 mm). A cSS band up to 2 voxels thick
    erodes away; a 3-voxel band leaves a long flat sliver that fails the aspect test. Count sub-objects that are
    microbleed-like: >= MB_MIN_VOX voxels across, compact (longest / shortest extent in mm <= MB_MAX_ASPECT), small
    (equivalent diameter of the re-dilated piece < MB_DIAM_MM) and mostly surrounded by parenchyma (>= MB_PARENCHYMA
    of a 1-voxel shell in labelled tissue). Microbleeds under ~3 mm erode to one voxel and are not counted."""
    vox = np.asarray(vox, float)
    s3 = np.ones((3, 3, 3), bool); s1 = ndi.generate_binary_structure(3, 1)
    if vox[2] > 1.5 * min(vox[0], vox[1]): s1[:, :, 0] = s1[:, :, 2] = False      # in-plane cross only
    er = ndi.binary_erosion(mask, structure=s1)
    lab, n = ndi.label(er, structure=s3)
    good = 0
    for i, sl in enumerate(ndi.find_objects(lab), 1):
        if sl is None: continue
        nv = np.array([s_.stop - s_.start for s_ in sl]); ext = nv * vox
        if nv.max() < MB_MIN_VOX or ext.max() / ext.min() > MB_MAX_ASPECT: continue
        piece = ndi.binary_dilation(lab == i, structure=s1) & mask
        diam = float((6 * piece.sum() * np.prod(vox) / np.pi) ** (1 / 3))
        if diam >= MB_DIAM_MM: continue
        shell = ndi.binary_dilation(piece, structure=s3) & ~mask
        if shell.any() and tissue[shell].mean() >= MB_PARENCHYMA: good += 1
    if good >= MB_MIN_N:
        return float(min(1.0, 0.5 + 0.1 * (good - MB_MIN_N))), good, f"{good} compact parenchymal blobs: microbleed cluster?"
    return 0.0, good, ""


class MimicRejector:
    """Per-scan state (filtered phase, scanner sign convention) + evaluate() per candidate."""

    def __init__(self, base, subj, swi, brain, vox, si, ves=None):
        self.vox = np.asarray(vox, float); self.si = si; self.aff = swi.affine; self.brain = brain
        self.hp = None; self.conv = 0; self.phase_note = "no phase image (data/ID_phase.nii.gz)"
        p = f"{base}/data/{subj}_phase.nii.gz"
        if not os.path.exists(p): return
        swi_p = f"{base}/data/{subj}_swi.nii"
        ph = nib.load(p)
        if ph.shape[:3] != swi.shape[:3] or not np.allclose(ph.affine, swi.affine, atol=1e-3):
            self.phase_note = "phase image is not on the SWI grid - phase rule off"; return
        if not fresh(p, swi_p):
            self.phase_note = "phase image is older than the SWI - phase rule off (re-import with --phase)"; return
        P, unit = phase_to_radians(np.asarray(ph.dataobj, np.float32))
        if P.ndim == 4: self.phase_note = "4-D phase - phase rule off"; return
        self.hp = homodyne_highpass(P, self.vox, si, HP_SIGMA_MM)
        self.hp[~brain] = 0
        # scanner sign convention: veins are paramagnetic -> fit the dipole pattern of the whole vessel map
        if ves is None or ves.sum() < 200:
            self.hp = None; self.phase_note = "vessel map too small to fit the phase sign convention - phase rule off"; return
        pred = highpass_linear(dipole_field(ves.astype(np.float32), self.vox, self.aff), self.vox, si)
        a, t, r = fit_amplitude(self.hp[ves], pred[ves])
        self.conv_t, self.conv_r = t, r
        if not np.isfinite(t) or abs(t) < CONV_MIN_T:
            self.hp = None
            self.phase_note = f"phase sign convention unclear on the vessel map (t {t:.1f}) - phase rule off"; return
        self.conv = 1 if a > 0 else -1
        self.phase_note = f"phase used ({unit}); paramagnetic sign {self.conv:+d} fitted on the vessel map (t {t:.1f}, r {r:.2f})"

    def phase_evidence(self, mm, pad):
        """-> (phase_z, phase_r) for the candidate mask mm (crop at slices pad); NaN without phase"""
        if self.hp is None: return np.nan, np.nan
        m = int(np.ceil(BOX_MARGIN_MM / self.vox.min()))
        idx = np.argwhere(mm) + [s.start for s in pad]
        lo = np.maximum(idx.min(0) - m, 0); hi = np.minimum(idx.max(0) + m + 1, self.hp.shape)
        box = tuple(slice(a, b) for a, b in zip(lo, hi))
        cm = np.zeros(tuple(hi - lo), bool); q = idx - lo; cm[tuple(q.T)] = True
        aff = self.aff.copy(); aff[:3, 3] = (self.aff @ np.r_[lo, 1])[:3]
        pred = highpass_linear(dipole_field(cm.astype(np.float32), self.vox, aff), self.vox, self.si)
        # fit on the shell AROUND the candidate: inside a strong deposit the phase wraps past pi (a paramagnetic sheet
        # flat to B0 at 0.5 ppm reaches 5.5 rad at 3 T, TE 20 ms) and the homodyne filter folds it to the wrong sign;
        # the field outside is weaker. Voxels already near +-pi are left out as possibly wrapped.
        reg = ndi.binary_dilation(cm, iterations=FIT_SHELL_VOX) & ~cm & self.brain[box] \
            & (np.abs(self.hp[box]) < WRAP_RAD)
        a, t, r = fit_amplitude(self.hp[box][reg], pred[reg])
        return (self.conv * t if np.isfinite(t) else np.nan), (self.conv * r if np.isfinite(r) else np.nan)

    def evaluate(self, mm, pad, tissue_pad, fz, fcz, dark_pad=None):
        """one candidate -> dict of the mimic columns. dark_pad: dark voxels (detector z < Z_LOW) in the same crop.
        The microbleed rule looks at the dark voxels within 2 voxels of the candidate, not at the candidate mask
        itself: the detector's mask is the ridge-filtered (thinned) dark region, in which round microbleeds no
        longer look round (PHM phantom: a row of 4 microbleeds became one candidate with 0 countable blobs)."""
        pz, pr = self.phase_evidence(mm, pad)
        support = mm if dark_pad is None else (mm | (ndi.binary_dilation(mm, iterations=2) & dark_pad))
        calc = 0.0; why = []
        if np.isfinite(pz) and pz <= CALC_Z and np.isfinite(pr) and pr <= CALC_R:
            calc = float(min(1.0, 0.7 + 0.05 * (CALC_Z - pz)))
            why.append(f"diamagnetic on phase (phase_z {pz:.1f}, r {pr:.2f}): calcification?")
        fs, fcls, fwhy = flair_rule(fz, fcz)
        ms, mn, mwhy = mb_cluster_rule(support, tissue_pad, self.vox)
        if fwhy: why.append(fwhy)
        if mwhy: why.append(mwhy)
        cls, sc = max((("calcification", calc), (fcls, fs), ("mb_cluster", ms)), key=lambda t: t[1])
        return dict(mimic_score=round(sc, 2), mimic_class=cls if sc > 0 else "none", mimic_reason="; ".join(why),
                    phase_z=round(pz, 2) if np.isfinite(pz) else np.nan,
                    phase_r=round(pr, 2) if np.isfinite(pr) else np.nan, mb_cluster_n=int(mn),
                    _flair_class=fcls)
