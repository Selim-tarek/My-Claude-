# v5.2 — surface-aligned patches, patch classifier, ranker upgrade

Source: the user's master plan (surface-aligned patches, `classify_surface_patches.py`, ablations A–G, required
tests) and Priority 3 of the v5.1 review prompt (nested CV, PR-AUC, Brier, baselines, model card).
Research software: everything here **ranks candidates for expert review**. Nothing is a diagnosis, nothing is a
calibrated probability, and nothing changes which objects are shown or excluded.

**Output-preserving.** The detector's candidates, exclusions, order and every existing CSV column are unchanged.
`detect_css.py` additionally writes `work/ID_pial_sd.nii.gz` (the signed pial distance it already computed), which the
patch extractor needs. With no model trained, nothing downstream changes.

## What was added

| # | Change | Where | Changes output? |
|---|---|---|---|
| 1 | **Surface-aligned patches.** Per object (shown and excluded) a patch of fixed physical size (default 24 × 12 × 8 mm, 0.5 mm grid) in a right-handed frame: u = the object's long axis projected into the cortical tangent plane, v = across the surface, n = outward pial normal (tissue → CSF) from the gradient of `work/ID_pial_sd.nii.gz`, mapped to world space with the inverse transpose of the affine (oblique / anisotropic data handled). The sign is checked against the labels and flipped if reversed (`normal_flipped`). Channels swi / phase / flair / labels on one world grid; a missing, off-grid or stale channel is **not zero-filled** — it is absent, its validity mask is all False and the reason is recorded. Intensity voxels outside the brain mask are never valid. Objects with < 3 voxels are listed as skipped with the reason (no silent loss). | `scripts/surface_aligned_patch.py` | new files `review/ID_patches/*.npz`, `review/ID_patches.csv` |
| 2 | **Patch classifier** (cSS vs vein): `features` (per channel in a fixed order: validity fraction, robust-z distribution, mean gradient along u / v / n, 8-bin profiles along each axis; label channel: tissue share along n; every channel adds `<ch>_missing`, `<ch>_valid_fraction`; optional geometry and vessel columns from the detector CSVs), `train` (dev split only, `--split test` refused, ≥ `--min-studies` studies and both classes; LOPO with nested L2; ablations A swi, B +phase, C +labels/geometry, D +vessel, E +flair, F swi+phase+labels+geometry+vessel, G all; AUC with patient-bootstrap CI, PR-AUC, Brier, sensitivity at 1 FP/scan; model JSON with feature schema, channel order, preprocessing, training studies, fold log, out-of-fold predictions, label hash, build, model hash, `research_only`), `predict` (refuses training studies, schema / channel-order mismatches and a model whose hash does not match; writes provenance per row). | `scripts/classify_surface_patches.py` | new files only |
| 3 | **Shared model code**: class-balanced L2 logistic regression, LOPO by patient with optional nested L2 selection (grid 0.1 / 1 / 10), one-class training folds skipped and logged (never filled in), AUC, average precision, Brier, patient bootstrap, sensitivity at k FP per scan, model hash. numpy / scipy only. | `scripts/css_model.py` | — |
| 4 | **Ranker upgrade**: `train_ranker.py` uses `css_model`; nested LOPO by default (`--no-nested` = fixed `--l2`); reports model, `score_v4`, darkness-only and vessel-only models on the same folds (AUC + CI, PR-AUC, sensitivity at 1 and 3 FP/scan, Brier for fitted models); writes `review/ranker_model_card.md`; the JSON adds `evaluation`, `lopo_fold_status`, `build`, `research_only`, `model_hash`. All keys `rank_candidates.py` reads are unchanged. | `scripts/train_ranker.py` | model JSON / card only |

## Tests (`tests/check_v52.py`, in `run_tests.sh`)

Master-plan tests covered: sheet vs tube patch (dark extent across the surface 12 mm vs 2 mm; sheet below, tube
above the surface; tissue on −n); reversed normal flipped by the labels; oblique 1 × 1 × 2 mm affine (n = rotated
slice axis); < 3 voxels skipped with reason; missing phase; FLAIR off the SWI grid and FLAIR older than the SWI
(each recorded as missing, features NaN, `<ch>_missing` 1); partial patch at the array edge (fixed shape, valid
fraction 0.68); outside-brain voxels never valid; patient grouping (a patient's out-of-fold scores ignore its own
labels, plain and nested); one-class fold skipped; schema persistence (predict refuses a schema mismatch); model hash
(a tampered model is refused); test split refused; too few studies refused; training study not predicted. On PH3 every
shown and excluded object gets a manifest row.

### External contract suite

`tests/verify_updated_architecture.py` (supplied by the user, unchanged, in `run_tests.sh`): 5 PASS, 1 SKIP. To meet
it, `surface_aligned_patch.extract_patch()` became the library entry (the CLI now calls it; PH3 patches are
bit-identical before / after) and `classify_surface_patches.py` gained `read_manifest` (duplicate or missing patch ->
error), `build_matrix` (features named `<channel>:<feature>`, e.g. `phase:missing`; labels other than cSS / vein codes
-> error), `grouped_lopo` and the `lopo-manifest` subcommand. SKIP: the optional `scripts/object_ledger.py` - the fate
invariant is already enforced inside `detect_css.py` (exit 4, manifest `fates`), so no separate module was added.

## What the synthetic test shows about the ablations (honest)

On the synthetic tables the single informative column sits among ~190 patch features: LOPO AUC A_swi 0.997, G_all
0.83 with 120 rows from 6 studies. **Adding channels costs accuracy when labelled studies are few** (more
weights, same data). The ablation table is there to show whether phase / FLAIR / geometry add anything on real
labels; do not default to G_all.

## Not done / limits

- No real labelled data yet: every number above is synthetic. The classifier must beat `score_v4` on held-out
  **test** studies (`evaluate.py`) before it is used for any reading order.
- The patch frame depends on the label-derived pial surface: on the SWI-only route (SynthSeg on SWI) it is coarse.
- Patch features are hand-designed summaries (no CNN); a learned encoder is for when ~10+ labelled patients exist.
- The 0.5 mm patch grid oversamples typical 0.8–1.2 mm SWI; patches carry the native resolution's blur.
