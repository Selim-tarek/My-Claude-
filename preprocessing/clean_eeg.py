"""
Resting-state EEG cleaning for brain-age prediction (ds005385 / Dortmund Vital Study).

Drop-in replacement for the preprocessing cells of `notebook_brainage_finale_luglio.ipynb`.

What this fixes relative to the notebook version
------------------------------------------------
1. Operation order. Average reference is applied *last*, after bad channels are
   detected and interpolated. Referencing before that smears one dead electrode
   across all 64 channels.
2. Bad channels are actually handled. `read_raw_bids` marks `status == "bad"` from
   channels.tsv into `raw.info["bads"]`; the notebook never interpolated them, so
   they flowed straight into the covariances, the PSDs and the deep-net window
   tensor. Here they are detected (BIDS + LOF + flat/noise heuristics), interpolated,
   and counted in QC.
3. Ocular / muscle artifacts are removed with ICA. Blink rate and blink amplitude
   both change with age, and blinks are far more frequent in eyes-open than
   eyes-closed. Leaving them in means the model can read age off eye movement
   rather than off brain activity - which would invalidate the EO/EC comparison
   in goal G3.
4. Real AutoReject. The notebook called `autoreject.get_rejection_threshold()`,
   which is only the global peak-to-peak threshold - not AutoReject, which repairs
   channels per epoch by interpolation instead of discarding the whole epoch.
5. High-pass moved to 0.5 Hz so the 1/f fit range no longer starts inside the FIR
   transition band (see `FOOOF_RANGE` note below). The aperiodic exponent and
   offset are the notebook's headline aging biomarkers, and a high-pass corner
   sitting at the low edge of the fit range biases both.
6. The epoch cache lives outside the run folder, so `FRESH_RUN = True` no longer
   silently discards hours of preprocessing (and re-downloads 74 GB).
7. Failures are recorded per recording instead of being swallowed by `except: pass`.

Notes on parameters that interact with the rest of the notebook
---------------------------------------------------------------
* `FOOOF_RANGE`: with `l_freq = 0.5` and MNE's default FIR design, the transition
  band spans roughly 0.25-0.75 Hz. Fitting the aperiodic component from 1.0 Hz is
  still close to that shoulder, so the default here is (2.0, 40.0). Changing the
  fit range changes the absolute exponent values, so keep it fixed across the
  whole cohort and report it.
* `epoch_len` stays at 10.0 s for drop-in compatibility with the deep nets
  (`n_times` is baked into ShallowFBCSPNet / the JEPA patch grid). A 10 s epoch
  plus a peak-to-peak rejection rule only works because ICA has already removed
  blinks; without ICA a 10 s window almost always contains one, and rejection
  becomes a filter on how much a subject blinks.

This module has no runtime dependency on the notebook's globals.
"""

from __future__ import annotations

import os
import json
import traceback
from dataclasses import dataclass, field, asdict
from typing import Iterable, Sequence

import numpy as np


# =====================================================================
# Configuration
# =====================================================================

@dataclass
class CleanConfig:
    # --- filtering ---
    sfreq: float = 250.0          # target sampling rate (downsample only, never up)
    l_freq: float = 0.5           # FIR high-pass edge (Hz)
    h_freq: float = 45.0          # FIR low-pass edge (Hz)
    notch: float | None = 50.0    # line noise (Germany); skipped when >= h_freq
    ica_l_freq: float = 1.0       # ICA is fit on a 1 Hz high-passed copy (Winkler 2015)

    # --- epoching ---
    epoch_len: float = 10.0
    epoch_overlap: float = 0.0

    # --- bad channels ---
    detect_bad_channels: bool = True
    lof_neighbors: int = 20
    flat_uv: float = 0.5          # channel std below this (uV) counts as flat
    noisy_z: float = 4.0          # robust z of log channel std above this counts as noisy
    max_bad_frac: float = 0.20    # more bads than this -> recording is unusable

    # --- ICA ---
    run_ica: bool = True
    ica_n_components: float = 0.99   # variance kept, or an int component count
    ica_method: str = "picard"       # falls back to "fastica" if python-picard missing
    ica_max_iter: int | str = "auto"
    eog_proxies: tuple[str, ...] = ("Fp1", "Fp2", "AF7", "AF8", "Fpz")
    ica_eog_threshold: float = 3.0
    ica_muscle_threshold: float = 0.8
    max_ica_reject: int = 6          # never drop more than this many components

    # --- muscle annotation (needs raw sfreq high enough for a 110-140 Hz band) ---
    annotate_muscle: bool = True
    muscle_z: float = 5.0
    muscle_min_good: float = 0.2

    # --- epoch rejection ---
    use_autoreject: bool = True
    ar_n_interpolate: tuple[int, ...] = (1, 4, 8)
    ar_consensus: tuple[float, ...] = (0.2, 0.4, 0.6, 0.8, 1.0)
    ar_fit_max_epochs: int = 120     # AutoReject is fit on a subsample, applied to all
    fixed_reject_uv: float = 200.0   # fallback peak-to-peak threshold

    # --- minimum usable data ---
    min_epochs: int = 6              # per condition, after rejection
    min_clean_seconds: float = 60.0  # per condition, after rejection

    # --- misc ---
    montage: str = "standard_1005"
    random_state: int = 42
    n_jobs: int = 1

    def as_dict(self) -> dict:
        return asdict(self)


#: Recommended specparam / FOOOF fit range given `CleanConfig.l_freq = 0.5`.
FOOOF_RANGE = (2.0, 40.0)


# =====================================================================
# Bad-channel detection
# =====================================================================

def detect_bad_channels(raw, cfg: CleanConfig) -> list[str]:
    """Union of BIDS-declared bads, flat channels, outlier-variance channels and LOF.

    Runs on filtered data. Returns channel names; does not modify `raw`.
    """
    import mne

    bads = set(raw.info["bads"])                    # from channels.tsv via mne-bids
    picks = mne.pick_types(raw.info, eeg=True, exclude=[])
    names = [raw.ch_names[i] for i in picks]
    data = raw.get_data(picks=picks)                # volts

    # flat / near-flat
    sd = data.std(axis=1)
    bads |= {n for n, s in zip(names, sd) if s < cfg.flat_uv * 1e-6}

    # variance outliers, via a robust z on log-std (immune to the flat channels above)
    finite = sd > 0
    if finite.sum() > 4:
        log_sd = np.log(sd[finite])
        med = np.median(log_sd)
        mad = np.median(np.abs(log_sd - med)) or 1e-12
        z = 0.6745 * (log_sd - med) / mad
        outliers = np.asarray(names)[finite][np.abs(z) > cfg.noisy_z]
        bads |= set(outliers.tolist())

    # local outlier factor on the channel-by-channel correlation structure
    try:
        lof = mne.preprocessing.find_bad_channels_lof(
            raw, n_neighbors=cfg.lof_neighbors, picks="eeg"
        )
        bads |= set(lof)
    except Exception:
        pass  # MNE < 1.7, or too few channels; the heuristics above still apply

    return sorted(bads)


# =====================================================================
# ICA
# =====================================================================

def _fit_ica(raw, cfg: CleanConfig):
    import mne

    raw_ica = raw.copy().filter(
        cfg.ica_l_freq, None, fir_design="firwin", verbose=False
    )
    method = cfg.ica_method
    try:
        import picard  # noqa: F401
    except Exception:
        if method == "picard":
            method = "fastica"

    ica = mne.preprocessing.ICA(
        n_components=cfg.ica_n_components,
        method=method,
        max_iter=cfg.ica_max_iter,
        random_state=cfg.random_state,
    )
    ica.fit(raw_ica, verbose=False)
    return ica, raw_ica


def _ica_exclude(ica, raw_ica, cfg: CleanConfig) -> tuple[list[int], dict]:
    """Pick ocular (and, when detectable, muscle) components. Returns (indices, detail)."""
    detail: dict = {"eog": [], "muscle": [], "eog_proxy": None}

    proxies = [c for c in cfg.eog_proxies if c in raw_ica.ch_names]
    if proxies:
        try:
            eog_idx, _ = ica.find_bads_eog(
                raw_ica, ch_name=proxies, threshold=cfg.ica_eog_threshold, verbose=False
            )
            detail["eog"] = sorted(set(int(i) for i in eog_idx))
            detail["eog_proxy"] = proxies
        except Exception:
            pass

    try:
        mus_idx, _ = ica.find_bads_muscle(
            raw_ica, threshold=cfg.ica_muscle_threshold, verbose=False
        )
        detail["muscle"] = sorted(set(int(i) for i in mus_idx))
    except Exception:
        pass  # requires MNE >= 1.1

    # Rank by explained variance so the cap keeps the most impactful components.
    combined = sorted(set(detail["eog"]) | set(detail["muscle"]))
    if len(combined) > cfg.max_ica_reject:
        try:
            var = ica.get_explained_variance_ratio(
                raw_ica, components=combined, ch_type="eeg"
            )["eeg"]
            var = np.atleast_1d(var)
            order = np.argsort(var)[::-1]
            combined = [combined[i] for i in order[: cfg.max_ica_reject]]
        except Exception:
            combined = combined[: cfg.max_ica_reject]
        combined = sorted(combined)

    detail["excluded"] = combined
    return combined, detail


# =====================================================================
# The pipeline for one recording
# =====================================================================

def clean_recording(bids_path, cfg: CleanConfig, qc: dict | None = None):
    """Read one BIDS recording and return cleaned, epoched, artifact-rejected `Epochs`.

    Order of operations:
        read -> pick EEG -> montage -> muscle annotation (wideband, pre-filter)
        -> notch -> band-pass -> downsample
        -> bad-channel detection -> interpolation -> AVERAGE REFERENCE
        -> ICA (ocular/muscle) -> fixed-length epochs -> AutoReject

    `qc` is filled in place with per-recording diagnostics. Returns `None` when the
    recording is unusable, with the reason recorded in `qc["reject_reason"]`.
    """
    import mne
    from mne_bids import read_raw_bids

    qc = qc if qc is not None else {}
    qc.setdefault("file", os.path.basename(str(getattr(bids_path, "fpath", bids_path))))

    raw = read_raw_bids(bids_path, verbose=False)
    raw.load_data()
    raw.pick("eeg")
    qc["n_channels_raw"] = len(raw.ch_names)
    qc["sfreq_raw"] = float(raw.info["sfreq"])
    qc["duration_raw_s"] = float(raw.times[-1]) if len(raw.times) else 0.0
    qc["bads_from_bids"] = list(raw.info["bads"])

    # --- montage: required for interpolation, RANSAC/LOF and every topoplot ---
    try:
        has_montage = raw.get_montage() is not None
    except Exception:
        has_montage = False
    if not has_montage:
        raw.set_montage(cfg.montage, on_missing="ignore")
    pos = raw.get_montage()
    if pos is None:
        qc["reject_reason"] = "no montage / no channel positions"
        return None
    missing_pos = [
        ch["ch_name"]
        for ch in raw.info["chs"]
        if not np.all(np.isfinite(ch["loc"][:3])) or np.allclose(ch["loc"][:3], 0)
    ]
    if missing_pos:
        # Channels without coordinates cannot be interpolated or plotted; drop them
        # rather than letting them silently poison the average reference.
        qc["dropped_no_position"] = missing_pos
        raw.drop_channels(missing_pos)

    # --- muscle annotation, before the 45 Hz low-pass destroys the 110-140 Hz band ---
    if cfg.annotate_muscle and raw.info["sfreq"] > 300:
        try:
            annot, _ = mne.preprocessing.annotate_muscle_zscore(
                raw, ch_type="eeg", threshold=cfg.muscle_z,
                min_length_good=cfg.muscle_min_good,
                filter_freq=(110, 140), verbose=False,
            )
            raw.set_annotations(raw.annotations + annot)
            qc["muscle_annot_s"] = float(np.sum(annot.duration)) if len(annot) else 0.0
        except Exception:
            qc["muscle_annot_s"] = None

    # --- spectral filtering ---
    if cfg.notch and cfg.notch < cfg.h_freq:
        # Only meaningful when line noise falls inside the pass-band; with h_freq=45
        # the 50 Hz notch is a no-op and is skipped.
        harmonics = np.arange(cfg.notch, min(cfg.h_freq, raw.info["sfreq"] / 2), cfg.notch)
        if len(harmonics):
            raw.notch_filter(harmonics, verbose=False)
    raw.filter(cfg.l_freq, cfg.h_freq, fir_design="firwin", verbose=False)

    # --- downsample (never upsample: that invents information) ---
    if raw.info["sfreq"] > cfg.sfreq + 1e-3:
        raw.resample(cfg.sfreq, verbose=False)
    elif raw.info["sfreq"] < cfg.sfreq - 1e-3:
        qc["reject_reason"] = f"sfreq {raw.info['sfreq']:.0f} Hz below target {cfg.sfreq:.0f} Hz"
        return None

    # --- bad channels: detect -> interpolate -> only THEN average reference ---
    if cfg.detect_bad_channels:
        raw.info["bads"] = detect_bad_channels(raw, cfg)
    qc["bads"] = list(raw.info["bads"])
    qc["n_bads"] = len(raw.info["bads"])
    if qc["n_bads"] > cfg.max_bad_frac * len(raw.ch_names):
        qc["reject_reason"] = f"{qc['n_bads']}/{len(raw.ch_names)} bad channels"
        return None
    if raw.info["bads"]:
        raw.interpolate_bads(reset_bads=True, verbose=False)

    raw.set_eeg_reference("average", projection=False, verbose=False)

    # --- ICA for ocular / muscle components ---
    qc["n_ica_excluded"] = 0
    if cfg.run_ica:
        try:
            ica, raw_ica = _fit_ica(raw, cfg)
            excl, detail = _ica_exclude(ica, raw_ica, cfg)
            ica.exclude = excl
            ica.apply(raw, verbose=False)
            qc["n_ica_components"] = int(ica.n_components_)
            qc["n_ica_excluded"] = len(excl)
            qc["ica_detail"] = detail
            del raw_ica
        except Exception as exc:
            qc["ica_error"] = f"{type(exc).__name__}: {exc}"

    # --- epoching ---
    epochs = mne.make_fixed_length_epochs(
        raw, duration=cfg.epoch_len, overlap=cfg.epoch_overlap,
        preload=True, reject_by_annotation=True, verbose=False,
    )
    qc["n_epochs_before_reject"] = len(epochs)
    if len(epochs) == 0:
        qc["reject_reason"] = "no epochs after annotation rejection"
        return None

    # --- epoch-level artifact handling ---
    epochs = _reject_epochs(epochs, cfg, qc)
    if epochs is None or len(epochs) == 0:
        qc.setdefault("reject_reason", "all epochs rejected")
        return None

    qc["n_epochs"] = len(epochs)
    qc["clean_seconds"] = len(epochs) * cfg.epoch_len
    qc["frac_epochs_kept"] = len(epochs) / max(1, qc["n_epochs_before_reject"])
    return epochs


def _reject_epochs(epochs, cfg: CleanConfig, qc: dict):
    """AutoReject (fit on a subsample, applied to all epochs), with a threshold fallback."""
    if cfg.use_autoreject:
        try:
            from autoreject import AutoReject

            ar = AutoReject(
                n_interpolate=np.array(cfg.ar_n_interpolate),
                consensus=np.array(cfg.ar_consensus),
                random_state=cfg.random_state,
                n_jobs=cfg.n_jobs,
                verbose=False,
            )
            n = len(epochs)
            fit_idx = (
                np.linspace(0, n - 1, cfg.ar_fit_max_epochs).astype(int)
                if n > cfg.ar_fit_max_epochs
                else np.arange(n)
            )
            ar.fit(epochs[np.unique(fit_idx)])
            cleaned, log = ar.transform(epochs, return_log=True)
            qc["reject_method"] = "autoreject"
            qc["n_epochs_interpolated"] = int((log.labels == 2).any(axis=1).sum())
            return cleaned
        except Exception as exc:
            qc["autoreject_error"] = f"{type(exc).__name__}: {exc}"

    epochs.drop_bad(reject=dict(eeg=cfg.fixed_reject_uv * 1e-6), verbose=False)
    qc["reject_method"] = f"fixed {cfg.fixed_reject_uv:.0f} uV p2p"
    return epochs


# =====================================================================
# Per-subject driver (pools runs within a condition)
# =====================================================================

def clean_subject_condition(bids_paths: Sequence, cfg: CleanConfig, qc_rows: list | None = None):
    """Clean every recording of one subject x condition and concatenate them.

    Pooling pre- and post-task runs maximises the data behind each covariance
    estimate while keeping all of a subject's data on one side of every CV split.
    """
    import mne

    parts = []
    for bp in bids_paths:
        rec_qc: dict = {}
        try:
            ep = clean_recording(bp, cfg, rec_qc)
        except Exception as exc:
            # Recorded, not swallowed: a silent `except: pass` turns a systematic
            # failure into an unexplained drop in subject count.
            rec_qc["error"] = f"{type(exc).__name__}: {exc}"
            rec_qc["traceback"] = traceback.format_exc(limit=3)
            ep = None
        if qc_rows is not None:
            qc_rows.append(rec_qc)
        if ep is not None:
            parts.append(ep)

    if not parts:
        return None
    if len(parts) == 1:
        return parts[0]

    try:
        parts = mne.equalize_channels(parts, verbose=False)
    except Exception:
        pass
    return mne.concatenate_epochs(parts, verbose=False)


def cache_path(cache_root: str, sub: str, cond: str) -> str:
    """Epoch cache path.

    Keep `cache_root` OUTSIDE the notebook's run folder. The notebook stored
    derivatives under RUN_DIR while `FRESH_RUN = True` renames RUN_DIR on every
    execution, so every run re-downloaded and re-preprocessed the whole cohort.
    """
    os.makedirs(cache_root, exist_ok=True)
    return os.path.join(cache_root, f"sub-{sub}_{cond}-epo.fif")


def config_fingerprint(cfg: CleanConfig) -> str:
    """Stable hash of the cleaning settings, written next to the cache.

    Guards against the classic silent bug: changing a filter setting and then
    reusing epochs cached under the previous settings.
    """
    import hashlib

    blob = json.dumps(cfg.as_dict(), sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:12]


def assert_cache_matches(cache_root: str, cfg: CleanConfig) -> None:
    """Raise if the cache on disk was written with different cleaning settings."""
    os.makedirs(cache_root, exist_ok=True)
    stamp = os.path.join(cache_root, "clean_config.json")
    fp = config_fingerprint(cfg)
    if os.path.exists(stamp):
        with open(stamp, encoding="utf-8") as fh:
            prev = json.load(fh)
        if prev.get("fingerprint") != fp:
            raise RuntimeError(
                f"Epoch cache in {cache_root} was written with a different CleanConfig "
                f"({prev.get('fingerprint')} != {fp}). Delete the cache or point "
                f"cache_root somewhere new."
            )
    else:
        with open(stamp, "w", encoding="utf-8") as fh:
            json.dump({"fingerprint": fp, "config": cfg.as_dict()}, fh, indent=2, default=str)


# =====================================================================
# Band-limited covariances (the input to the Riemannian models)
# =====================================================================

def band_covariances(epochs, bands: dict, sfreq: float, picks: Sequence[str] | None = None):
    """Per-band OAS covariance, filtering each epoch independently.

    The notebook concatenated all epochs into one long array and filtered that,
    which puts a filter transient at every 10 s epoch boundary - the joins are not
    continuous recordings. Filtering per epoch (with reflection padding) and only
    then concatenating for the shrinkage fit removes those transients.
    """
    import mne
    from sklearn.covariance import OAS

    ep = epochs.copy().pick(list(picks)) if picks is not None else epochs
    data = ep.get_data()                                   # (n_epochs, C, T)
    n_ch = data.shape[1]
    covs = np.empty((len(bands), n_ch, n_ch))

    for bi, (lo, hi) in enumerate(bands.values()):
        xb = mne.filter.filter_data(
            data, sfreq, lo, hi, pad="reflect_limited", verbose=False
        )                                                   # filtered per epoch
        flat = np.concatenate(list(xb), axis=1)             # (C, n_epochs*T)
        covs[bi] = OAS().fit(flat.T).covariance_
    return covs


# =====================================================================
# QC audit: is data quality itself an age proxy?
# =====================================================================

def qc_age_audit(qc_df, age_col: str = "age", verbose: bool = True):
    """Correlate every QC quantity with age.

    This is the check the notebook is missing. If retained-epoch count, bad-channel
    count or the number of removed ICA components correlates with age, then the
    model can reach a low MAE by reading data quality rather than brain state, and
    the drivers reported for G3/G4 are contaminated. Older subjects genuinely tend
    to produce noisier recordings, so a non-zero correlation is expected - the point
    is to measure and report it, and to consider matching or regressing it out.
    """
    import pandas as pd
    from scipy import stats

    cols = [
        c for c in qc_df.columns
        if c != age_col and pd.api.types.is_numeric_dtype(qc_df[c]) and qc_df[c].notna().sum() > 3
    ]
    rows = []
    for c in cols:
        d = qc_df[[age_col, c]].dropna()
        if d[c].nunique() < 2 or len(d) < 4:
            continue
        r, p = stats.pearsonr(d[age_col], d[c])
        rho, p_s = stats.spearmanr(d[age_col], d[c])
        rows.append(dict(metric=c, n=len(d), pearson_r=r, pearson_p=p,
                         spearman_rho=rho, spearman_p=p_s))

    out = pd.DataFrame(rows)
    if len(out):
        out = out.reindex(out["pearson_r"].abs().sort_values(ascending=False).index)
        out = out.reset_index(drop=True)
    if verbose:
        print("QC-vs-age audit (|r| descending):")
        print(out.to_string(index=False) if len(out) else "  (no numeric QC columns)")
        if len(out) and out["pearson_p"].min() < 0.05:
            worst = out.iloc[0]
            print(
                f"\n  [!] '{worst['metric']}' correlates with age "
                f"(r={worst['pearson_r']:+.2f}, p={worst['pearson_p']:.3g}). "
                f"Report this, and check whether the brain-age gap survives "
                f"controlling for it."
            )
    return out


__all__ = [
    "CleanConfig",
    "FOOOF_RANGE",
    "detect_bad_channels",
    "clean_recording",
    "clean_subject_condition",
    "cache_path",
    "config_fingerprint",
    "assert_cache_matches",
    "band_covariances",
    "qc_age_audit",
]
