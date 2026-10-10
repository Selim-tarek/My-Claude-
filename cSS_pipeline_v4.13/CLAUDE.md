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
raw/<ID>/SWI, raw/<ID>/T1, raw/<ID>/FLAIR   dcm2niix output (series-numbered NIfTI + JSON)
data/<ID>_swi.nii           skull-stripped SWI used by everything downstream (never the minIP)
data/<ID>_t1.nii.gz, data/<ID>_flair.nii.gz   imported T1 / FLAIR (prep_anat.sh)
synthseg/<ID>_seg.nii.gz    SynthSeg --parc --robust on the SWI (fallback when no T1)
work/<ID>_t1seg_swispace.nii.gz   T1 labels (SynthSeg or recon-all aparc+aseg) on the SWI grid
work/<ID>_a2009s_swispace.nii.gz  Destrieux labels on the SWI grid (recon-all only) -> sulcal scoring
work/<ID>_flair_swispace.nii.gz   FLAIR on the SWI grid
work/<ID>_swi2t1.lta, <ID>_flair2swi.lta   registrations
work/<ID>_seg_swispace.nii.gz  labels used by detect (T1-derived if present, else SWI SynthSeg)
work/<ID>_dark.nii.gz       dark-voxel map used for region growing at scoring
work/<ID>_ich.nii.gz        optional reader-drawn lobar ICH; work/<ID>_ich_used.nii.gz = mask used
review/<ID>_candidates.{csv,nii.gz}  ranked candidates (cand_id = rank)
review/<ID>_calls.csv       reader decisions + candidate fingerprint (vox_i/j/k, volume)
review/css_scores.csv       one row per subject
subjects/<ID>/              FreeSurfer recon-all (SUBJECTS_DIR)
logs/<ID>_run.log           run_all.sh output
**Never commit patient data. IDs are study codes (P006…), never names/MRNs.**

## Pipeline
run_css.sh ID swi.nii [--t1 T1] [--flair FLAIR] [--recon] → SynthStrip (CSF kept) →
  prep_anat.sh (T1 SynthSeg or recon-all + bbregister, FLAIR coreg) or SynthSeg on SWI →
  align_seg.py → detect_css.py → freeview
review_css.py ID (all candidates) → reader calls → mark_css.py → score_css.py
Validation: make_synthetic.py → eval_synthetic.py → stress_test.py; feature_report.py on real calls.

## detect_css.py (v4) — key logic (full description: PIPELINE.md §7)
- pial surface: signed distance to labelled tissue (− in tissue, + in CSF), normals from its gradient
- zone = pial band (−3.5..+4 mm) minus WM, deep nuclei and a 3 mm skull-strip rim
- darkness z = (I − local cortex mean, σ 10 mm, 2nd pass excludes cortex z<−2) / cortical MAD
- per-slice Sato ridge; strong z<−2.5 & ridge>85th pct; weak z<−2.0 & >80th; hysteresis per hemisphere
- filters: ≥25 mm³, ≥2 slices, elongation ≥2.0 (tram-track ≈2.4)
- definition features: pial_dist_mm, bank_frac, tube_ratio (3-D Hessian), surface_alignment,
  tram_frac, vein_tree_mm (recorded only: merges brain-wide on real SWI), mirror_dark_frac, ich_dist_mm/near_ich_suggest, infratentorial,
  flair_csf_z/flair_bright
- score = score_v4 (untrained definition rule) unless CSS_RANK=v3; score_v3 always recorded
- CSV columns up to long_structure are consumed by export_review.py and the workbook — keep stable;
  new columns are appended after them.

## score_css.py (v4)
- Destrieux available: foci → sulci (≤4 mm, ≥15 % of voxels); sulci adjacent if touching or
  bordering a common gyrus; per hemisphere 0 / 1 (≤3 adjacent sulci) / 2; STRIVE-2 by sulci
- otherwise Euclidean approximation (3 mm foci, 10 mm clusters)
- infratentorial and reader near-ICH excluded from 0–4 but reported; n_unreviewed reported
- ICH rule (Charidimou 2017, verified from the PDF): ≥3 unaffected sulci between cSS and lobar ICH
  (≥2 at multiple axial levels if no superficial path) — counted on the Destrieux adjacency graph;
  <2 excluded with a reader-drawn ICH mask (warned with the auto mask), ==2 flagged
- sequence (SWI / T2*-GRE from raw/<ID>/SWI/*.json ImageType) and % sulci affected per hemisphere
  (van Harten 2023) are recorded
- Boston v2.0 cSS count by gyri (Lancet Neurol 2022): gyri adjacent if touching or sharing a sulcus;
  1 component of <=2 gyri = 1 lesion, else >=2 (boston2_css_lesions)
- detect: cmb_like = extent <=10 mm and >=50 % of the outer shell in parenchyma (AJNR 2016 rule)
- optional SWI phase (run_css.sh --phase -> data/<ID>_phase.nii.gz): reviewer shows homodyne
  high-pass phase (calcium = opposite sign to veins); score records TE, field, vendor, voxel size

## Validation status (honest)
- v4 first real run: P006 with T1 SynthSeg (no recon-all yet): 71 candidates; top ranks in R inferior
  parietal and L supramarginal (the areas of the v3-accepted cSS); not yet reviewed. vein_tree_mm
  was 255 mm for most candidates (network merges brain-wide) -> made recorded-only.
  Auto-ICH flagged skull-base artifact + sagittal sinus (reader: no ICH, alignment good) -> ICH
  blobs must now be >3 mm under the pial surface, off-midline, off the rim, outside artifact zones.
- P006 full route (recon-all, 100 candidates): detector found 6/7 of the v3-accepted cSS (ranks 5-49;
  missed one R precentral = iron-rich motor cortex?). Same reader, blinded re-review: accepted 1/100
  (called the old cSS Normal/Vein) -> non-expert labels unreliable; expert labels needed.
  Scoring the 7 old foci: L=1 (fixed vs v3 L=2, expert L1), R=1 (expert R2): R sulci intraparietal,
  postcentral, superior temporal were linked via a shared (supramarginal) gyrus -> "shared gyrus"
  adjacency may be too loose; --explain added to inspect links before changing the rule (n=1).
- Phantoms: v4 phantom (tram-track/convexity/single-bank cSS vs tubular sulcal and cortical
  surface veins, ICH, Destrieux labels): 6/6 found, every definition feature separates cSS from
  veins, ICH Dice 0.86, sulcal score 3/4 as constructed. v3 phantom: 8/8 synthetic lesions; its
  "veins" are flat planes in sulcal banks (cSS-shaped), so its ranking AUC is not meaningful for v4.
- Synthetic (UCSF CMB_labeler test scans P001–P005, v3): strong ~92–100 %, medium ~83 %, faint ~46–63 %.
- Real case P006 (Mayo, 1.2 mm SWI, expert R2+L1=3), v3: pipeline+reader 4 (L2 vs 1 from the
  Euclidean adjacency rule). feature_report (v3, 17 calls): score AUC 0.53, darkness 0.37 (veins
  darker) → v4 damps darkness (√). Reader was not blinded for P006.

- P006 --explain (reader's 7 old foci): R: IPS-postcentral touch (7 mm); IPS-STS via angular gyrus
  (16 mm); postcentral-STS via planum temporale (21 mm) -> one group, R=1. L: subcentral-STS via
  supramarginal (27 mm) -> L=1. A distance cut cannot give L1 and R2 (27 mm adjacent vs 16 mm not);
  the v3 R precentral focus (not detected in v4) would make R 4 sulci -> R2 = expert. Rule NOT
  changed; waiting for expert labels (expert_sheet.py / expert_import.py).
- v4.7: css_evidence / vein_evidence (rule-based 0-1 summaries, not probabilities); mirror recorded
  only; score records TR, slice thickness, phase availability.

- v4.8 FP reduction: merge pieces <3 mm apart; exclude specks <6 mm, microbleed-like, tubular
  candidates >0.5 mm out in sulcal CSF without tram-track, skull-base artifact zone (lowest 30 %),
  edge contrast <0.75 SD. Excluded -> review/ID_excluded.csv/.nii.gz; CSS_KEEP_ALL=1 disables.
  Phantoms: 6/6 and 8/8 kept, PH3 veins shown 6 -> 3, PH1S candidates 14 -> 12 (AUC .69 -> .75).
  score_css sulcus assignment: >=15 % of the focus OR >=20 mm3 (merged foci keep both sulci).
  check_known.py reports whether known lesions survive a detector change.

- v4.9 normal anatomy: exclude candidates with >=30 % within 5 mm of cerebellum/brainstem (tentorial
  interface), >=30 % within 5 mm of brainstem/ventral DC/hippocampus/amygdala (basal cisterns), or in
  the lowest 20 % of the cerebrum height (artifact labels: lowest 40 %). Columns infra_frac,
  basal_frac, rel_height. PH4 test phantom (cerebellum + dark surface line): line excluded, 6/6 kept.
  Trade-off: true cSS on the inferior temporal/occipital surface over the tentorium is excluded too
  (listed in ID_excluded.csv).

- v4.9.1: merging capped at 40 mm (P006: unlimited 3 mm merging chained veins into 10-17 cm objects,
  one swallowed known focus #13 and was excluded as tentorial). PH1S synthetic now 7/8: one low
  synthetic lesion falls in the skull-base band (lowest 20 % of the cerebrum) - known trade-off.

- P010 (cSS-negative SWI-only partial slab, ~0.47 mm in-plane): 67 FP after 33 excluded; stress test
  sensitivity strong 56 % / medium 25 % / faint 19 %, ranking AUC ~0.5, features ~0.5 (SynthSeg-on-SWI
  anatomy too coarse). check_known on P010S: lesions excluded as basal cisterns / skull base /
  microbleed-like. v4.9.2: skull-base + basal rules need T1 anatomy AND a scan that includes the
  skull base (coverage check); SWI-only infratentorial needs >=50 %; cmb_like flag only; merge gap
  2x slice (3-4 mm); elongation filter only for single-piece candidates (merged V-shapes were lost).
  Phantoms back to 8/8, 16/16, 6/6. Conclusion: T1 should be part of the study protocol.

- v4.9.2 on P010 (SWI-only): strong 12/16 (75 %), medium 4/16, faint 3/16; clean FPs 99 (position
  rules off); ranking/features ~0.5. P006 (T1): 70 candidates, 6/6 known kept (ranks 8-45).
  Conclusion: hand rules at their limit; SWI-only geometry features carry no signal.
- v4.9.3: stress_test.py uses the host's T1 labels when present (copies work/H_t1seg_swispace to
  HS after the synthetic SWI is written). Next: stress test on P006 = do features work WITH T1?

- P006 stress (T1, v2 synthetic): strong 94 %, medium 94 %, faint 50 %; ranking AUC .87/.78/.61;
  features bank_frac .85, branch .78, surface_gradient .73, contact .72, volume .74, elongation .38
  -> partly artifacts of the v2 generator (big, blob-like lesions).
- v4.10 make_synthetic v3: thin traced lines (1-2 vox, 5-30 mm log-uniform, 2-5 slices), tram-track
  40 %, vein decoys (--veins, in sulcal CSF) -> work/HS_veins.nii.gz; stress_test --veins 4 default,
  --legacy; feature table adds AUC vs decoy veins.

- P011 (cSS-negative, T1 SynthSeg on a 2.5 mm reformat): 26 candidates (26 excluded); stress v3 strong 94 %,
  medium 69 %, faint 31 %, AUC ~0.95; features bank_frac .86, cortex_dist .21, mirror .12, vein_tree .24.
- v4.11 vessel exclusion: vessel map = dark Hessian tubes in the cortical band + WM up to 15 mm (>=3 SD
  darker than WM, pieces >=3 mm, round blobs and ICH rim excluded, 1-voxel gaps bridged). Features
  vessel_wm_mm (tube continuing straight on into deep WM inside a 30 deg cone from the candidate's
  deep end), axis_normal (|long axis . normal|), vessel_run_mm (straight-on continuation beyond the ends,
  recorded only), vessel_ext_mm (recorded only). Rule "vessel": wm >=4 mm AND axis >=0.5 AND tram <0.4.
  Pitfalls found on the phantom: axis alone removes cSS wrapping crown->sulcus; un-coned WM reach
  connects to noise / ICH rim. Phantom: 24/24 synthetic kept, transcortical veins excluded.
  make_synthetic --radial K (transcortical decoys), stress_test reports missed-lesion reasons and decoy fate.
  Not yet checked on real data: run check_known on P006, stress on P011.

- v4.11 on real data: P006 6/6 known kept, 70 -> 66 candidates (4 "vessel"); P011 unchanged (26, no
  vessel exclusions); P011 stress 94/75/25 %, AUC .94/.98/.96, missed lesions 18 not detected + 3 "vein",
  0 "vessel" (rule is safe but rarely fires). Decoys mostly not detected (31/36) -> decoy test
  uninformative on P011. Real FPs lie ALONG the surface: axis_normal higher in synthetic cSS (AUC .82,
  synthetic closed-sulcus artifact), vessel_wm .46, vessel_run .32. vein_tree (.23) and mirror (.12)
  separate synthetic from real FPs, but known P006 cSS also have vein_tree 255 -> partly synthetic artifact.
- rule_test.py: per-feature thresholds just beyond the most extreme KNOWN lesion (P006 6 foci) -> FPs
  removed on negatives (P011), synthetic lesions lost. Shortlist only (n=6).

- rule_test P006 (6 known) vs P011 (26 FP): best safe rules score_v4<1.03 (22 FP, 4/34 synth lost),
  pial_dist>-0.08 / cortex_dist>0.27 / cortex_frac<0.46 (20 FP, 6-8/34 synth), bank_frac<0.45 (14, 2/34),
  edge<2.17 (14, 3/34). Known cSS: cortex_frac >=0.55, pial <=-0.57, mirror median .17 (FP .93), all
  score_v4 >=1.24. v4.12 adopts a CONSERVATIVE definition rule "off the cortex": T1 only, cortex_frac<0.30
  AND pial_dist>0.30 AND tram<0.2. score threshold not adopted (ranking, would cut synthetic).

- v4.12 on real data: P006 66 -> 59 candidates (9 "off the cortex"), 6/6 known kept; P011 26 -> 8
  (18 "off the cortex"); P011 stress unchanged 94/75/25 %, AUC .93/.98/.97, 0 synthetic lesions lost to the
  new rule. Remaining P011 FPs: 5 orbitofrontal (frontal-sinus susceptibility zone, rel_height >= 0.4),
  rostral ACC, precuneus, middle temporal. Candidate next rule: orbitofrontal artifact zone - needs more
  negatives first.

- v4.13 from van Harten 2023 (full text re-read incl. Discussion): score_css adds surface_mm2 and
  L/R_surface_pct (cortex voxels facing CSF within 1.5 mm of the grown mask; voxel face area = V^(2/3))
  + volume_note; phantom blooming test +1 vox: volume x2.13 vs surface x1.30. agreement.py: Pearson,
  ICC(A,1), Bland-Altman, Dice between two review-folder copies; --vs-score (their Fig. 4).
  Not done: follow-up/progression mode (register follow-up SWI to baseline, new cSS volume/area/sulci).

- v4.14 Stage 1 of docs/AUDIT_v4.13_and_FP_plan.md (bookkeeping, no change to shown/excluded candidates on
  the phantoms - old CSV columns identical in values on PH1/PH1S/PH2/PH3/PH4/PH5): css_qc.py -> review/ID_qc.json
  (protocol, orientation, coverage, anatomy, registration cost from prep_anat work/ID_reg.json, label-pial
  vs SWI-edge agreement overall + L/R x inferior/middle/superior, background ratio) and status ok /
  low_confidence_review / insufficient_quality - PROVISIONAL thresholds, REPORTED ONLY. Appended columns
  qc_status, rule_version (E1-E9 ids), n_pieces; excluded.csv adds exclusion_confidence. Generation-gate
  drops (volume/slices/elongation/no cortex) now logged in work/ID_dropped.csv/.nii.gz; stress_test names
  them in "missed synthetic lesions". Fix B1: extent_mm was 2-3.5 mm too small for objects on array slice 0
  (partial slabs -> false "speck" exclusions). Phantom: edge agreement 0.12 / 0.88 / 1.38 / 1.62 mm for
  label shifts 0 / 0.8 / 1.6 / 3.2 mm (only the 3.2 mm shift crosses the 1.5 mm threshold). Not yet run on
  P006/P011.

- v4.15 expert-label tooling (plan section 4; no detector change): css_common.LABELS = fine codes C, V/VS/VC/VT/VX
  (surface / sulcal / transcortical / sinus veins), TV, A/AM, N/NI, MB, CA, SAH, LN, HI, IS, H, U, each mapped to
  the old coarse call (cSS/Vein/Normal/Artifact/Near ICH/Unsure) so review, feature_report, mark/score are
  unchanged; optional :1-:5 confidence. expert_sheet --include-excluded mixes excluded objects in, unmarked.
  expert_import --session N writes labels/candidate_labels.csv (long table incl. excluded objects, defaulted
  flag) and warns if an excluded object is called cSS. label_missed.py: reader-drawn outlines of cSS with no
  number -> labels/missed_lesions.csv with fate shown / excluded (rule id) / dropped (gate) / not generated.
  study_meta.py: labels/study_meta.csv, patient-level reference, QC/protocol fields from ID_qc.json, split
  dev/test fixed once assigned (--force logged). P006/P010/P011 must stay 'dev'. tests/check_labels.py.

## Known issues / decisions
- SWI-only route (SynthSeg on SWI) gives a coarse pial surface → pial features less reliable.
- v4 ranking weights are hand-set from the definition, not trained; review all candidates.
- Medial cSS is mildly down-ranked by the midline factor (0.7).
- QSM/SEPIA not used (site decision). Phase only as a reading aid (calcium mimic); it cannot separate
  hemosiderin from deoxy-Hb.
- minIP series must not be used as input; the reviewer computes its own minIP slab for reading.
- Workbook (tools/build_workbook.py) builds the Excel validation workbook; recalculated with LibreOffice.

## Next steps
1. Re-run P006 with --t1 --recon (back up review/P006_* first), blinded re-review; target L1 R2 = 3.
2. feature_report after each case; keep features with consistent AUC ≥0.75 / ≤0.25.
3. ~15 cases: logistic regression on the definition features, leave-one-patient-out CV.
4. Second blinded reader (κ), cSS-negative controls for false-positive burden.

## How to test
`bash tests/run_tests.sh` — v3 phantom (radial sulci, fissure, veins) end to end, fails if synthetic
sensitivity < 4/8; v4 phantom (check_v4.py: sensitivity, score_v4 AUC ≥0.8, ICH Dice, sulcal 3/4).
tests/check_labels.py (v4.15 label tools) and tests/check_qc.py (both called by run_tests.sh): v4.14 QC status, column contract, dropped log, extent at the
volume edge, label-shift monotonicity, transposed-orientation detection.
`bash tests/test_shell.sh` — run_css.sh / prep_anat.sh / run_all.sh flow with stubbed FreeSurfer.
No FreeSurfer or patient data needed.

## Conventions
- Python 3.11, numpy/scipy/nibabel only (+matplotlib for the reviewer); no network calls.
- Every script: CLI with `sys.argv`, reads `$CSS_BASE`, prints a one-line summary first.
- Prefer small, verifiable changes; run tests/run_tests.sh after any change to detect/score logic.
