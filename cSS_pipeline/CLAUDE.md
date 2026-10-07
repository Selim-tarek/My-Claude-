# CLAUDE.md — context for Claude Code

## Project
Semi-automatic grading of **cortical superficial siderosis (cSS)** on SWI MRI, for cerebral amyloid
angiopathy research (Mayo Clinic, Jacksonville; cerebrovascular/neuroimaging fellowship project).
Output: candidate cSS foci for human review → validated **cSS multifocality score 0–4**
(per hemisphere: 0 none; 1 = one sulcus or ≤3 immediately adjacent sulci; 2 = ≥2 non-adjacent or >3
adjacent sulci; total = left + right) plus continuous measures (foci, volume, regions).
Research software only — not for clinical use.

## Environment
macOS (Apple M3, 36 GB), FreeSurfer 8.2.0 (`mri_synthseg`, `mri_synthstrip`, `recon-all`, `freeview`),
Miniforge conda env `css` (python 3.11, nibabel, numpy, scipy, scikit-image, pandas, matplotlib).
All scripts read the data root from `$CSS_BASE` (default `~/css_project`).

## Data layout ($CSS_BASE)
raw/<ID>/SWI, raw/<ID>/T1   dcm2niix output (series-numbered NIfTI + JSON)
data/<ID>_swi.nii           skull-stripped SWI used by everything downstream
synthseg/<ID>_seg.nii.gz    SynthSeg --parc --robust on the SWI (1 mm grid)
work/<ID>_seg_swispace.nii.gz  labels resampled to SWI grid, masked to SWI coverage
work/<ID>_dark.nii.gz       dark-voxel map used for region growing at scoring
review/<ID>_candidates.{csv,nii.gz}  ranked candidates (cand_id = rank)
review/<ID>_calls.csv       reader decisions from review_css.py
review/css_scores.csv       one row per subject
subjects/<ID>/              FreeSurfer recon-all (SUBJECTS_DIR)
**Never commit patient data. IDs are study codes (P006…), never names/MRNs.**

## Pipeline
run_css.sh ID swi.nii → SynthStrip (CSF kept) → SynthSeg → align_seg.py → detect_css.py → freeview
review_css.py ID → reader calls → mark_css.py → score_css.py
Validation: make_synthetic.py → eval_synthetic.py → stress_test.py; feature_report.py on real calls.

## detect_css.py (v3.2) — key logic
- zone = (cortex dilated 2 in-plane | 2 mm outer edge) minus 3 mm 3-D rim (skull-strip edge artifact)
- darkness z = (I − local cortex mean, Gaussian σ 10 mm) / global robust SD (MAD), after 99.5 pct clip
  (v3.2: LOCAL reference — global reference flagged normal iron-rich motor cortex bilaterally)
- per-slice Sato ridge (σ 1–3 vox, black ridges); strong: z<−2.5 & ridge>85th pct; weak: z<−2.0 & >80th
- hysteresis: weak components kept only if they contain a strong voxel; labelled per hemisphere
- filters: ≥25 mm³ (regionprops with spacing already returns mm³ — earlier versions double-scaled),
  ≥2 slices, elongation ≥2.5
- score = darkness × (0.25+surface_contact) × log1p(vol); ×0.5 artifact zone / midline <5 mm /
  long structure >40 mm in z; ×0.6 vein_like (cortex_frac<0.35)
- recorded only: cortex_dist_mm, branch_per10mm, surface_gradient
- CSV column names are consumed by export_review.py and the Excel workbook — keep them stable.

## Validation status (honest)
- Synthetic (UCSF CMB_labeler test scans P001–P005: 3T GE, skull-stripped 2 mm slabs, radiation-induced
  microbleeds, presumed cSS-negative): strong lesions ~92–100 %, medium ~83 %, faint ~46–63 % found.
  Synthetic lesions are drawn on the cortical boundary, so some features are partly circular there.
- Real case P006 (Mayo, 1.2 mm SWI, expert score R2+L1=3): cSS present in top 20 (first at rank 3);
  pipeline+reader score 4 (R2 ✓, L2 vs 1). Left discrepancy = two accepted foci (inferior parietal,
  supramarginal) 36 mm apart → Euclidean 10 mm adjacency rule calls them non-adjacent.
- feature_report on P006 (17 calls: 7 cSS, 5 vein, 4 normal, 1 artifact): ranking score AUC 0.53;
  darkness 0.37; cortex_frac 0.43; surface_gradient 0.24 (opposite to hypothesis); n_slices 0.73.
  → current features do NOT separate cSS from veins on real data. n=1: do not re-tune weights yet.
- Reader was not blinded to the expert score for P006.

## Known issues / decisions
- SynthSeg-on-SWI cortex is too coarse for vein-vs-cSS geometry → move anatomy to T1/FreeSurfer.
- Scoring adjacency is Euclidean (3 mm merge, 10 mm adjacency) — approximation of "adjacent sulci".
- Near-ICH cSS is excluded from 0–4 (reader flag `--ich`) but count/volume are reported.
- Medial cSS is down-ranked by the midline penalty (known trade-off).
- QSM/SEPIA not used (site decision). Phase images exist but are not used.
- Workbook (tools/build_workbook.py) builds the Excel validation workbook; recalculated with LibreOffice.

## Next steps (v4 — needs recon-all)
1. Register SWI→T1: `bbregister --s ID --mov data/ID_swi.nii --reg work/ID_swi2t1.lta --t2 --init-coreg`;
   QC overlay in freeview.
2. Bring candidates to T1 (`mri_vol2vol --nearest`), compute per candidate:
   signed distance to pial surface (veins in CSF above pial; cSS on it), mirror check
   (contralateral homologous darkness via surface registration or midsagittal flip), vein-network
   connectivity to the superior sagittal sinus, thickness variation along skeleton.
3. Sulcal scoring: assign foci to Destrieux (aparc.a2009s) sulcal labels; two sulci adjacent if they
   border the same gyrus on the pial mesh; replace Euclidean grouping in score_css.py.
4. Re-run P006 (target: left = 1, total 3) and feature_report; then learned classifier
   (e.g. logistic regression on features + reader calls) after ~15 reviewed cases.

## How to test
`bash tests/run_tests.sh` — phantom brain (radial sulci, fissure, veins), runs every script end to end,
fails if synthetic sensitivity < 4/8. No FreeSurfer or patient data needed.

## Conventions
- Python 3.11, numpy/scipy/nibabel only (+matplotlib for the reviewer); no network calls.
- Every script: CLI with `sys.argv`, reads `$CSS_BASE`, prints a one-line summary first.
- Prefer small, verifiable changes; run tests/run_tests.sh after any change to detect/score logic.
