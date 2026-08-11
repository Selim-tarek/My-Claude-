# EEG cleaning review — `notebook_brainage_finale_luglio.ipynb`

Review of the preprocessing stage (cells 6–16) of the ds005385 brain-age notebook,
plus a corrected implementation in [`clean_eeg.py`](clean_eeg.py).

The modelling half of the notebook is careful — subject-level folds, de Lange & Cole
correction fit on training data only, tangent-space reference fit inside CV. The
cleaning half is where the problems are, and cleaning errors propagate into every
model downstream, so a bad pipeline can produce a good-looking MAE for the wrong
reason.

---

## Findings, most consequential first

### 1. Bad channels are never repaired

`read_raw_bids` copies the `status` column of `channels.tsv` into
`raw.info["bads"]`. The notebook's `preprocess_recording` never reads that field,
never detects bads itself, and never interpolates. A dead or drifting electrode
therefore goes straight into:

- the OAS covariances (the input to every Riemannian model),
- the Welch PSD (and so the `specparam` 1/f exponent and offset),
- the deep-net window tensor.

`set_eeg_reference("average")` does exclude `info["bads"]` when computing the
reference — but the bad channel itself stays in the data, and every subsequent
stage consumes it as if it were signal.

**Fix:** detect (BIDS status ∪ flat ∪ variance outliers ∪ LOF) → interpolate →
then reference.

### 2. Average reference is applied before filtering and before bad-channel repair

```python
raw.set_eeg_reference("average", ...)   # ← first
raw.notch_filter(...)
raw.filter(1.0, 45.0, ...)
```

Averaging first spreads drift, line noise, and whatever a bad electrode is doing
across all 64 channels. Once it is in the reference it cannot be removed from the
other channels. Correct order is filter → repair bad channels → average reference.

### 3. No ocular artifact removal, and this one is not cosmetic

There is no ICA anywhere in the notebook. For this specific study that is a
validity problem, not just a noise problem:

- Blink rate and blink amplitude both change with age.
- Blinks are far more frequent in eyes-open than eyes-closed.

Goal **G3** asks which scalp locations and which task (EO/EC) drive the age
prediction. With blinks left in, a frontal/prefrontal EO effect is exactly what
you would expect from ocular artifact alone — the analysis cannot distinguish
"prefrontal EEG ages" from "older subjects blink differently". Any frontal or
EO-specific finding is unfalsifiable until blinks are removed.

**Fix:** ICA fit on a 1 Hz high-passed copy, ocular components identified against
Fp1/Fp2 as EOG proxies (this dataset has no dedicated EOG channel), muscle
components via `find_bads_muscle`, capped at 6 removed components.

### 4. `get_rejection_threshold()` is not AutoReject

```python
thr = autoreject.get_rejection_threshold(epochs, decim=2)
epochs.drop_bad(reject=thr)
```

The markdown promises "AutoReject's data-driven" rejection. This function only
computes a global peak-to-peak threshold by cross-validation. Actual `AutoReject`
learns per-channel thresholds and **repairs** epochs by interpolating the few worst
channels instead of discarding the whole epoch — which matters a lot at 10 s
epochs, where one bad channel costs you the entire window across all 64 channels.

### 5. 10 s epochs + a 200 µV peak-to-peak rule is a bad combination

At a typical resting blink rate, most 10 s eyes-open windows contain at least one
blink, and a blink is 100–200 µV. So with the fallback threshold the notebook
either rejects most EO data or (if you raise the threshold) rejects nothing
meaningful. Worse, what survives is a function of *how much the subject blinked* —
and blink rate varies with age, so rejection itself becomes an age-correlated
filter.

The combination only becomes safe once ICA has removed the blinks first, which is
how the replacement pipeline orders it. Epoch length stays at 10 s so the deep nets
(`ShallowFBCSPNet`'s `n_times`, the JEPA patch grid `T // PT`) do not need retuning.

### 6. The 1 Hz high-pass sits at the low edge of the 1/f fit range

`CONFIG["l_freq"] = 1.0` and `CONFIG["fooof_range"] = (1.0, 40.0)`. With MNE's
default FIR design a 1 Hz high-pass has a transition band of roughly 0.5–1.5 Hz, so
the lowest points of the aperiodic fit sit **inside the filter roll-off**. The fit
is anchored on attenuated data, which biases both the exponent and the offset —
the two headline aging biomarkers in the notebook's own framing (§3, G3, G5).

**Fix:** high-pass at 0.5 Hz, fit the aperiodic model from 2.0 Hz
(`clean_eeg.FOOOF_RANGE`). Absolute exponent values shift with the fit range, so
fix it once for the whole cohort and report it.

### 7. The 50 Hz notch is dead code

`notch = 50.0` with `h_freq = 45.0`: the low-pass removes the line component
before the notch could matter, and the 45 Hz low-pass applied before `resample(250)`
already prevents aliasing. Harmless, but it is not doing what the markdown says it
does. The replacement only applies a notch when the line frequency actually falls
inside the pass-band.

### 8. `FRESH_RUN = True` destroys the epoch cache every run

```python
CONFIG["deriv_root"] = os.path.join(RUN_DIR, "derivatives")   # cache lives here
if FRESH_RUN and os.path.isdir(RUN_DIR): os.rename(RUN_DIR, RUN_DIR + ".bak_...")
```

The cached epochs, covariances, PSDs and `specparam` fits all live under `RUN_DIR`,
which is renamed away at the top of every execution. So every "Run All" re-streams
74 GB from S3 and redoes the whole preprocessing — on a run the notebook itself
estimates at 10–30 h. Move the cache outside the run folder
(`clean_eeg.cache_path`), and stamp it with a config fingerprint
(`assert_cache_matches`) so you can never silently reuse epochs cleaned under
different settings.

### 9. Silent failures

```python
try:
    e = preprocess_recording(bp)
    ...
except Exception:
    pass
```

A systematic failure — one channel naming convention, one corrupt file type —
shows up only as a smaller subject count with no explanation. Every failure is now
recorded per recording with its exception type and a short traceback.

### 10. Covariance filtering crosses epoch boundaries

```python
X = np.concatenate(list(data), axis=1)      # (C, n_epochs*T)
xb = mne.filter.filter_data(X, sfreq, lo, hi)
```

Epochs are not contiguous in time (rejected ones are gone), so every join is a
discontinuity, and filtering across it injects a transient — worst in delta, where
the filter is longest. `band_covariances()` filters each epoch independently with
reflection padding, then concatenates for the shrinkage fit.

### 11. No QC-versus-age audit

The notebook writes `qc_epoch_counts.csv` and never analyses it. If retained-epoch
count, bad-channel count, or removed-ICA-component count correlates with age, the
model can score well by reading data quality instead of brain state. Older cohorts
genuinely do produce noisier recordings, so some correlation is expected — the
point is to measure it and report it. `qc_age_audit()` does this and flags the
strongest correlate.

### Smaller points

- `if abs(sfreq - 250) > 1e-3: raw.resample(250)` will **upsample** a recording
  below 250 Hz, inventing bandwidth. The replacement downsamples only and rejects
  under-sampled recordings.
- `set_montage(..., on_missing="ignore")` leaves channels with no coordinates.
  Those cannot be interpolated or plotted, and they still enter the average
  reference. They are now dropped explicitly and logged.
- The common-channel intersection is computed by reading every epochs file twice,
  and one subject with an odd channel set silently shrinks the channel space for
  the entire cohort with no warning printed.

---

## Using the replacement

### Replacing cell 10

```python
import sys; sys.path.insert(0, "/path/to/repo")
from preprocessing.clean_eeg import (
    CleanConfig, FOOOF_RANGE, clean_subject_condition,
    cache_path, assert_cache_matches, band_covariances, qc_age_audit,
)

CLEAN = CleanConfig(
    sfreq=CONFIG["sfreq"],
    l_freq=0.5,                       # was 1.0 - keeps the 1/f fit out of the roll-off
    h_freq=CONFIG["h_freq"],
    epoch_len=CONFIG["epoch_len"],    # 10.0, unchanged, so the deep nets still fit
    run_ica=True,
    use_autoreject=RUN["autoreject"],
    fixed_reject_uv=CONFIG["fixed_reject_uv"],
    random_state=SEED,
)
CONFIG["fooof_range"] = FOOOF_RANGE   # (2.0, 40.0)

# Cache OUTSIDE RUN_DIR so FRESH_RUN does not throw away hours of work.
CACHE_ROOT = os.path.expanduser("~/eeg_cache/ds005385")
assert_cache_matches(CACHE_ROOT, CLEAN)
```

### Replacing cell 11

```python
qc, rec_qc = [], []
for sub in tqdm(subjects, desc="clean subjects"):
    rec = {"subject": sub, "age": age_of[sub]}
    cached = all(os.path.exists(cache_path(CACHE_ROOT, sub, c)) for c in ("EO", "EC"))
    if not cached:
        try:
            ensure_subject(sub)
        except Exception as e:
            print(f"  [download fail] sub-{sub}: {e}")
            rec.update(n_EO=0, n_EC=0); qc.append(rec); continue

    for cond, tlist in [("EO", eo_tasks), ("EC", ec_tasks)]:
        fp = cache_path(CACHE_ROOT, sub, cond)
        if os.path.exists(fp):
            rec[f"n_{cond}"] = len(mne.read_epochs(fp, verbose=False))
            continue
        rows = []
        ep = clean_subject_condition(get_recording_paths(sub, tlist), CLEAN, rows)
        for r in rows:
            r["subject"] = sub; r["condition"] = cond
        rec_qc += rows
        if ep is None or len(ep) < CLEAN.min_epochs:
            rec[f"n_{cond}"] = 0
            continue
        ep.save(fp, overwrite=True, verbose=False)
        rec[f"n_{cond}"] = len(ep)
        rec[f"bads_{cond}"]  = float(np.mean([r.get("n_bads", 0) for r in rows]))
        rec[f"ica_{cond}"]   = float(np.mean([r.get("n_ica_excluded", 0) for r in rows]))
        rec[f"kept_{cond}"]  = float(np.mean([r.get("frac_epochs_kept", np.nan) for r in rows]))
    qc.append(rec)
    if DELETE_RAW and not cached:
        shutil.rmtree(os.path.join(bids_root, f"sub-{sub}"), ignore_errors=True)

qc = pd.DataFrame(qc)
qc["has_both"] = (qc.get("n_EO", 0) > 0) & (qc.get("n_EC", 0) > 0)
pd.DataFrame(rec_qc).to_csv(os.path.join(CONFIG["deriv_root"], "qc_per_recording.csv"), index=False)
qc.to_csv(os.path.join(CONFIG["deriv_root"], "qc_epoch_counts.csv"), index=False)

# The check the original pipeline never ran:
qc_age_audit(qc).to_csv(os.path.join(CONFIG["deriv_root"], "qc_age_audit.csv"), index=False)
```

### Replacing `band_covariances` in cell 13

```python
def band_covariances(sub, cond):
    ep = mne.read_epochs(cache_path(CACHE_ROOT, sub, cond), verbose=False)
    from preprocessing.clean_eeg import band_covariances as _bc
    return _bc(ep, CONFIG["bands"], CONFIG["sfreq"], picks=common)
```

---

## Cost and expected effect

ICA plus AutoReject is the expensive part: roughly 1–3 min per recording on CPU,
so about 4–8× the current preprocessing time — but it is paid once, because the
cache now survives across runs. Set `run_ica=False, use_autoreject=False` for a
smoke test; the code path is otherwise identical.

Expect the reported MAE to get **worse**, not better, once blinks and bad channels
are removed. That is the point. A model reading blink rate and electrode drift will
score well on chronological age and mean nothing as a biomarker — a brain-age gap
is only interpretable if the signal it comes from is brain signal.

## Not verified at runtime

This container has no `numpy`/`mne`/`autoreject` installed and no access to the
ds005385 data, so `clean_eeg.py` is syntax-checked only — it has not been executed
against real recordings. The MNE calls that vary across versions
(`find_bad_channels_lof` needs ≥ 1.7, `ICA.find_bads_muscle` needs ≥ 1.1,
`get_explained_variance_ratio` needs ≥ 1.2) are each wrapped so an older version
degrades instead of failing. Run it on 2–3 subjects and check
`qc_per_recording.csv` before committing to the full cohort.
