# Changelog v5.1.0 (first stage of the v5.1 review)

Research software: expert review of every candidate remains required; the pipeline is not validated for clinical use.

Scope of this stage: the v5.1 review prompt's Priority 0 and the Priority 1 items that do not change which candidates
are shown, plus one bug fix in an opt-in v5.0.0 rule. With no switch set (`CSS_USE_MIMIC`, `CSS_QC_STRICT` unset),
candidates, exclusions and their order are unchanged against v5.0.0 / v4.19.3. The final test log
(`docs/TEST_LOG_v5.1.0.txt`) lists what ran.

## Changes

| # | Change | Files | Risk to sensitivity | Migration |
|---|---|---|---|---|
| 1 | **Fate invariant.** Every merged object ends in exactly one of dropped / excluded / shown, voxel for voxel; otherwise the detector prints the violation and exits **4** | detect_css.py | none (bookkeeping) | new exit code 4 |
| 2 | **Fate counts in the manifest**: raw components, not generated (weak component without a strong seed), pieces before merging, objects, dropped per gate, excluded per rule, shown, QC-blocked | detect_css.py | none | `review/ID_manifest.json` gains `detect.fates` |
| 3 | **Dropped log**: `gate_id` (G1_volume, G2_slices, G3_elongation, G4_no_cortex), `measured`, `threshold` appended | detect_css.py | none | new columns at the end of `work/ID_dropped.csv` |
| 4 | **Manifest schema** `css-manifest/1.0` and a **build id** (commit + date, filled in by `git archive` into `scripts/BUILD`; `git rev-parse` in a source checkout) | css_common.py, scripts/BUILD, .gitattributes, install.sh | none | `install.sh` copies BUILD |
| 5 | **Obliquity**: `obliquity_deg` (slice normal vs scanner S-I) and the array-to-world axis codes in the QC file; > 20 deg -> `low_confidence_review` ("oblique slices") | css_qc.py | none (QC status only) | oblique scans now read low confidence; with `CSS_QC_STRICT=1` that also scales `score_model` x0.8 |
| 6 | **Stale derived inputs**: surface features older than the detection are ignored by `train_ranker.py` / `rank_candidates.py` (`CSS_QC_STRICT=1`: stop). The model records the label-file hash and the detector builds of its training outputs; `rank_candidates.py` warns when the labels changed or the study was detected by another build (strict: refuses) | css_common.py (`stale_against`, `use_derived`), train_ranker.py, rank_candidates.py | none | models fitted before v5.1 have no hash: no check |
| 7 | **`audit_case.py SUBJECT`**: inputs (hash, geometry, orientation, obliquity, S-I axis), acquisition, freshness of T1 labels / FLAIR / phase / Destrieux / surface features, registration and QC, the fate counts and their consistency with the files -> printed + `review/ID_audit.json` | audit_case.py (new) | none | read-only |
| 8 | **evaluate.py**: per-study `obj_generated`, `obj_dropped`, `obj_excluded`, `obj_shown`, `obj_qc_blocked` from the manifest | evaluate.py | none | new columns at the end of `_patients.csv` |
| 9 | Repetition time in the acquisition record | css_qc.py | none | |
| 10 | **Bug fix, microbleed-cluster rule** (opt-in E11 of v5.0.0): with `CSS_USE_MIMIC=1` it **excluded true cSS lesion 3 of PH3** (a tram-track): eroding the dark support split the tram-track into compact fragments counted as microbleeds. Now separate blobs are counted without erosion (connected components that are compact, 2-6 mm, >= 2 voxels across in two directions, mostly surrounded by parenchyma) | mimic_rules.py | lower risk than v5.0.0 | `mb_cluster_n` values change (recorded column); touching microbleeds that form a chain are no longer counted |

## Already in place before v5.1 (verified by the new tests, not changed)

| v5.1 prompt item | Where |
|---|---|
| CSF = brain voxels that are not tissue (recon-all label 0 = SynthSeg label 24) | detect_css.py, v4.17; now tested (check_v51 part 2) |
| Orientation: the S-I axis is read from the affine; non-axial storage stops with `insufficient_quality` (exit 2) unless `CSS_ALLOW_ORIENT=1` | detect_css.py, v4.18; axis-swap invariance and oblique affines tested (check_v51 parts 3-4) |
| Sato / Hessian scales, merge gap, extent in millimetres | v4.18 (Sato mm), v4.9 / v4.14 |
| Hessian on the mask-filled image | v4.17 (audit B3) |
| `extent_mm` correct at the volume edge | v4.14 (audit B1, check_qc part 5) |
| Mirror samples ignore points outside the brain | v4.18 (`brain_in`) |
| Stale FLAIR / T1 labels / phase / Destrieux ignored | v4.17, v4.18.3, v5.0, score_css.py |

## Not done in this stage (planned)

- Per-piece features and `mixed_object_review` for merged objects (prompt P1.6).
- Valid-sample fractions for every reference statistic (P1.3, partly).
- Versioned expert-label schema with required confidence / adjudication status (P2.1).
- The feature families of `FEATURE_SPEC_cSS_vs_mimics.md` that do not exist yet (sulcal centre-line offset, compactness / sphericity, tortuosity, sinus connection, missingness indicators ...).
- Review-priority tiers (P2.3), the nested-CV ranker with PR-AUC / Brier / model card (P3), the surface-aligned patch extraction and patch classifier (master plan), richer phantoms (P4).

## Tests added

`tests/check_v51.py`: fate invariant on every phantom case; CSF label 0 vs 24 gives identical candidates and
exclusions; in-plane axis swap gives the same objects at the same world positions; 15 deg oblique affine changes no
candidate and 25 deg reads low confidence; stale FLAIR / phase / surface features / model labels; and **the mimic
layer switched on loses no true lesion on PH3 or PH6** (the check that would have caught item 10).

## v5.1.1 — rest of Priority 1 (per-piece evidence, valid-sample fractions, missing inputs, label schema)

Output-preserving with the defaults (every v5.1.0 column identical on the phantom suite; new columns appended).

| # | Change | Where | Default output |
|---|---|---|---|
| 1 | **Per-piece evidence of merged objects**: each pre-merge piece is classified with the detector's own vein_like definition (off the cortex > 0.3 mm and tubular > 0.4 vs on the cortex and not tubular); columns `piece_vol_vox`, `piece_pial_mm`, `piece_tube`, `piece_vein_frac`, `mixed_object` (>= 2 pieces, each kind >= 15 % of the volume) | css_pieces.py, detect_css.py | recorded only |
| 2 | **`CSS_MIXED_REVIEW=1`** (opt-in): a mixed object that a vein rule (E5 / E6 / E7) would exclude is kept for review, `mixed_object_review` = 1 | detect_css.py | off |
| 3 | **Valid-sample fractions**: `ref_valid_frac` (share of the local reference with usable cortex), `mirror_valid_frac` (mirror points inside the brain; below 0.5 -> `mirror_dark_frac` NaN, not "dark"), `region_n_vox` | detect_css.py | recorded only |
| 4 | **Missing-input indicators** `missing_t1_anatomy`, `missing_flair_features`, `missing_phase_features` (detector) and `missing_surface_features` (ranker); the ranker uses them as features BEFORE median imputation | detect_css.py, train_ranker.py, rank_candidates.py | recorded only |
| 5 | **Versioned expert labels**: `label_schema_version` (css-labels/1.0), `label_class`, `confidence_missing`; `expert_import.py --require-confidence` refuses unrated rows; consensus writes `adjudication_status` | css_common.py, expert_import.py, consensus.py | new label columns |

Limits: a vein TOUCHING a cSS segment voxel-to-voxel is already one piece (cannot be split); merges undone by the 40 mm
chain cap leave single pieces. A mixed phantom could therefore not be staged through the detector: the piece logic is
unit-tested (band + tube -> mixed, vein fraction 0.40; one kind / 12 % vein / no merging -> not mixed). PH3: 0 mixed
objects, `CSS_MIXED_REVIEW=1` changes nothing. Lopsided head mask: mirror_valid_frac 0.00-0.92 on 5 objects.
Bug found in the v5.1.1 full-suite run and fixed before commit: `add_missing` crashed when a study had no surface /
phase / FLAIR columns at all (`True.astype`).

## v5.1.2 — feature-spec families (FEATURE_SPEC_cSS_vs_mimics.md), recorded only

New module `css_features.py` (feature schema css-features/1.0), columns appended to both CSVs, no rule or score uses them:
- position: `pial_sd_p10/p90`, `surface_contact_frac_1mm/2mm`, `csf_adjacent_frac`, `csf_depth_mm` (distance into the
  CSF from tissue), `sulcal_depth_mm` (depth below a 5 mm closing of the tissue = gyral-crown envelope);
- compactness: `esd`, `compactness_3d` (marching-cubes surface), `sphericity`, `elongation_3d`, `flatness_3d`;
- path: `path_length_mm`, `tortuosity` (two-sweep Dijkstra on the largest 26-connected component);
- darkness: `darkness_z_p10`, `dark_frac_z2/3/4`; FLAIR `flair_csf_bright_frac`; phase `phase_valid_frac`,
  `phase_polarity_frac` (mimic_rules.py); `feature_qc_flags`, `feature_schema_version`.
Tests (check_v512): compactness sphere 1.0 / sheet 0.005 / tube 0.093; tortuosity straight 1.01, half circle 1.59;
PH3 directions of the spec all met (VT vessel_wm > cSS; VC / VS deeper in CSF than cSS; VS nearer the crowns than VC;
VT more compact than cSS); PHM phase valid fraction 1.00, polarity 0.57-0.96. Phantom directions only.

## v5.1.3 — review-priority tiers (reading order only)

`css_tiers.py`: `review_priority_rank` 1-5, `review_priority`, `review_priority_reason` appended; `review_css.py
--by-priority` (opt-in) reads in that order. 1 high-priority cSS candidate; 2 cSS-vs-mimic uncertain (mixed object,
possible calcification, FLAIR context); 3 likely vein — only on CONVERGENT evidence (WM continuation with a long dark
run or a tube; sinus zone with >= 10 mm continuation; tubular and continuing at both ends; off the cortex AND in the
CSF AND vessel-shaped), never one feature; 4 microbleed-like / air-bone artifact zone; 5 case QC insufficient. Every
object stays in its list. Cut-offs are placeholders (v4 rules + phantoms). PH3 / PH6: every cSS lesion tier 1, every
constructed vein tier 3; 16 unit rows (a single vein-ish feature alone stays tier 1). The tier-3 phantom criterion was
restated once during development (shown phantom veins had tube ratio <= 0.44 and first read tier 1): it is fitted to
phantoms and must be checked against expert labels (evaluate.py) before any use.

## v5.2

See docs/CHANGELOG_v5.2.md (surface-aligned patches, patch classifier, shared model code, ranker upgrade).
