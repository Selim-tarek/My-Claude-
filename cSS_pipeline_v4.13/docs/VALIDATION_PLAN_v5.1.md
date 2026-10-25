# Validation plan v5.1 (draft for the PI)

Research software: this plan describes how the pipeline would be evaluated. No result in it exists yet. As of v5.1.0
there are **no expert-labelled patients**; every number in the repository comes from phantoms or from one
development patient (P006) with an earlier, non-consensus read.

## 1. Frozen split

- Each study gets its split (`dev` or `test`) with `study_meta.py --split` **at enrolment, before anyone looks at the
  detector output**. The split cannot be changed without `--force`, which is logged.
- Development studies so far: P006, P010, P011 (rules were set on them) and the LEMON volunteers (P201, P202, P204 ...).
  They are never used to report performance.
- `train_ranker.py` refuses `--split test`; `run_css.sh` ranks only dev studies automatically; a test study is scored
  once, with the frozen model and frozen rules, then `evaluate.py --split test`.

## 2. Labels (reference standard)

- Fine codes (`css_common.LABELS`): C (cSS), VS / VC / VT / VX (surface / mid-sulcal / transcortical / sinus-adjacent
  vein), TV, MB, CA, SAH, LN / HI, A / AM, N / NI (normal / iron-rich cortex), H (ICH-related), IS, U (uncertain);
  each mapped to a coarse analysis class (`label_class`).
- Two blinded readers (`expert_sheet.py`, randomised order, no detector score shown), confidence 1-5, consensus with
  adjudication (`consensus.py`). Agreement is reported (kappa, `agreement.py`).
- The reference includes shown candidates, **excluded objects** (sheets with `--include-excluded`) and
  **reader-outlined missed cSS** (`label_missed.py`), so losses at generation, gates and exclusion rules are counted.
- To do (v5.1 plan): version the label schema and make confidence and adjudication status required fields.

## 3. Metrics (all patient-grouped, with patient-bootstrap 95 % intervals)

Primary:
- lesion-level and patient-level cSS sensitivity, overall and among the top 5 / 10 / 20 shown candidates;
- false-positive candidates per **cSS-negative** scan (reader-confirmed negatives only);
- agreement of the 0-4 multifocality score with the expert score (weighted kappa).

Secondary:
- FROC (rank cut-off -> sensitivity vs FP per scan); review burden (candidates to read per patient);
- reference cSS lost at each stage: not generated / dropped by gate / excluded by rule / QC-blocked (from the
  manifest fates and the missed-lesion outlines);
- for a learned ranker: PR-AUC, sensitivity at fixed FP per scan, calibration (Brier), compared on the same folds
  with score_v4, darkness only and vessel features only. Pooled candidate-level AUC is reported only as secondary.

Strata: sequence (SWI / GRE / SWI built from GRE), field strength, TE, in-plane resolution, slice thickness
(<= 1.2 / <= 2 / > 2 mm), anatomy route (recon-all / SynthSeg / SWI-only), QC status.

## 4. Sample-size rationale (to be decided with the PI)

Precision, not power: with patient-level sensitivity near 0.8, the 95 % interval half-width is about
1.96 x sqrt(0.8 x 0.2 / n): n = 27 cSS-positive patients gives +-0.15, n = 50 gives +-0.11. False positives per
negative scan need enough negatives for a stable mean: about 20-30 confirmed cSS-negative scans.
Hard-negative targets for developing mimic features (feature spec section 9): >= 25 VS, 25 VC, 15 VT, 15 CA, 25 MB
adjudicated objects. These are development targets, not a claim of adequate clinical sample size.

## 5. Acceptance criteria (pre-specified)

- A v5.x layer (mimic exclusion, ranker, tier) becomes a default only if, on held-out patients, it loses no more
  reference cSS than an agreed limit (proposal: no lesion-level sensitivity loss beyond the lower bound of the
  current interval) and reduces FP per negative scan with an interval that excludes zero change.
- No hard exclusion from a feature that has only phantom support (E10 / E11 stay opt-in until then).
- Every reported number is patient-grouped, stratified as above, and comes from studies whose split was frozen
  before the evaluated version existed.
- External validity is not claimed until a second site / vendor cohort has been evaluated.

## 6. Current status and next steps

1. Run the real-data regression checklist (RUNBOOK 13a) on P006 and P201 for the v5.x switches.
2. Expert reading of the development patients and confirmation of LEMON negatives.
3. Freeze the test split for newly enrolled patients and the label schema before any test-set reading.
