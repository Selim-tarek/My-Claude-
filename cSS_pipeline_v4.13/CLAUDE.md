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

- v4.16 surface-following veins (hard-negative class VS), all RECORDED ONLY - no rule, score or tier uses them:
  css_vessel.py local vein tracking in detect (cand_geodesic_mm, local_cont_max/min/total_mm, terminal_ends,
  local_net_mm, local_branch_points/_per10mm): geodesic from each candidate end through vessel or dark (>=3 mm
  pieces) voxels, forward 60 deg cone, <=20 mm, pial band only (WM noise faked continuation; noise voxels
  chained through gap bridging until dark pieces <3 mm were dropped). surface_features.py (recon-all only):
  candidates + excluded projected on lh/rh.pial; transform from orig.mgz + swi2t1.lta with the direction
  CHOSEN by agreement with the label pial boundary (stops if > 2 mm); surf_normal_dist_mm, surf_pieces,
  surf_band_continuity, geo_length_mm, geo_euclid_ratio, surf_dark_cont_mm, surf_exit_paths. Replaces the
  uploaded geodesic_vein_features.py prototype (its mesh "continuation" always hit the 25 mm limit and its
  "branches" were mesh vertices of valence >= 6). PH6 phantom (PH3 + branching surface vein -> sagittal sinus):
  vein continues >= 19 mm, cSS <= 7 mm; normal distance cSS -0.74..-0.85 mm vs vein/sinus +0.53..+1.43 mm;
  surf_pieces/dark_cont/exit_paths do NOT separate on the phantom. make_synthetic --surface-veins K
  (default 2 in stress_test, decoy fate per kind); finite random-walk decoys cannot test continuation.
  expert_sheet: coronal + sagittal panels; --vessel-overlay optional (off: detector output biases readers).
  surface_overlay.py (replaces the uploaded plot_css_pysurfer.py idea: no PySurfer/Mayavi, no unvalidated p_css
  or tiers): every object voxel -> nearest pial vertex (css_surface.py transform + 2 mm check), FreeSurfer
  annotation coloured by detector list or expert calls (safer class wins on overlap: cSS > uncertain > veins >
  artifact > other), id overlay (.mgh, +cand/-excl), per-object Destrieux sulcus/gyrus from
  label/?h.aparc.a2009s.annot, freeview command printed. tests/check_overlay.py.

- v4.17 audit Stage 2 + review follow-ups. Detector: CSF = brain & ~tissue (B2; confirmed on a phantom with
  recon-all labels: v4.13 surface_contact 0.00 for every candidate, FLAIR CSF reference = basal label-24 region),
  3-D Hessian on If (B3; phantom tube_ratio change <= 0.08, 0.17 on sinus pieces at the mask edge; same objects,
  same exclusions, ranks shift slightly), E8 needs T1 (B6), stale FLAIR / Destrieux / phase / LTA ignored or
  refused (B11), work/ID_gen.nii.gz generation-stage bits, review/ID_manifest.json (sha256 of inputs, params,
  package + FreeSurfer versions; surface_features adds its section). label_missed splits "not generated" into
  outside zone (rim / WM label / beyond band), not dark, no ridge, no strong seed. consensus.py: two-reader
  consensus, adjudication_todo.csv + --apply, Cohen kappa (cSS vs not, 6 classes) with patient-bootstrap CIs;
  agreement.py adds quadratic weighted kappa for multifocality / L / R score. environment.yml. Tests:
  check_stage2.py, check_consensus.py. Still to confirm on P006 recon-all labels: count of label 0 vs 24 in
  brain & ~tissue (the fix is safe either way).

- v4.17 on real data (user's Mac, FreeSurfer 8.2): P006 60 candidates, 6/6 known kept (ranks 6-43); P011 8
  candidates; P006 recon-all CSF outside tissue = label 0: 377,534 vox, label 24: 615,819 vox (B2 fix matters).
  surface_features crashed: real bbregister LTAs are type 0 with "cras" (parser expected "c_ras").
- v4.17.1: LTA parser reads "cras"; phantom test with a bbregister-style type-0 LTA. QC: P006/P011 edge
  agreement 1.62 / 2.12 mm mixed scatter with misregistration (P006 signed median -0.12, P011 +1.62). New
  iterative robust rigid fit (search window re-centred on the estimate, outliers > 1.5 mm dropped):
  label_shift_mm, label_shift_vec_mm, edge_bias_mm, edge_scatter_mm, edge_found_frac; status uses shift and
  bias only (> 1.5 mm). Phantom: true 0 / 1.13 / 2.56 / 5.12 mm -> 0.02 / 1.11 / 2.56 / 4.3 mm (the old fixed
  window read 5.12 mm as 1.28 mm and passed it).

- v4.17.1 on real data: P006 QC ok (shift 0.21 mm), surface LTA direction chosen 0.58 vs 15.51 mm, overlay 105
  objects. P011 low_confidence_review: label_shift_vec_mm [1.34, -0.07, -3.87] (found_frac 0.97) = 3.9 mm
  up-down slip from mri_coreg (no recon-all); axial freeview looked aligned. Re-registration with bbregister pending.
- v4.17.2: prep_anat.sh stops at the first failed step (a placeholder --t1 path fell through and kept the old
  labels); run_css.sh stops if prep_anat fails; --t1 data/ID_t1.nii.gz (the imported copy) works. RUNBOOK:
  coronal/sagittal registration check. test_shell 6a/6a2.
- v4.18 (external review, verified against the code): Sato sigmas 0.6/1.2/1.8 mm -> in-plane px (was 1-3 px:
  P010 0.47 mm vs P011 0.86 mm got different filters; values reproduce P006 at 0.573 mm); CSS_SATO_LEGACY=1 =
  old [1,2,3] px for before/after runs. mirror_dark_frac ignores mirror points outside the brain (NaN if most
  are). branch_per10mm uses true skeleton step length (css_common.skeleton_length). detect stops (exit 2, QC
  written) on a non-axial array, CSS_ALLOW_ORIENT=1 overrides. Phantoms: same candidate/excluded counts except
  PH1S (0.5 mm in-plane, sigmas become 1.2/2.4/3.6 px): 8/8 -> 7/8, lesion 3 (11 mm3 true) generated but its
  component measured 23.5 mm3 < 25 mm3 volume gate; stress small 16/16 -> 15/16 (same lesion); legacy switch
  restores 8/8 exactly. Fine-resolution scans (P010-like) are where the change matters most -> real-data gate:
  P006 check_known 6/6, P010 + P011 counts and stress test, both settings. docs/METHODS_scope_and_protocol.md
  drafts the phenotype-scope / protocol text (scope option a/b is a PI decision).
- v4.18 recorded-only additions: review_tier / tier_reason (2 = artifact zone, vein-like, cmb-like or vein_evidence >
  css_evidence; case QC not used), review_css --by-tier opt-in; css_vessel end_type_max/_min = stops / runs_on /
  branches / midline (continuation reaches dmid < 4 mm from an end that is not there). check_known prints both.
  Phantoms only: PH3/PH6 true lesions all tier 1, all non-lesions tier 2; PH6 surface veins run_on (19.8-20 mm),
  none typed branches/midline (the 20 mm reach stops before the sinus); unit geometry recovers all four types.
  Untested on real data: cSS next to veins may read runs_on/branches (cf. vein_tree 255 on known P006 cSS).
- v4.18 radius (recorded only; idea from Avadiappan 2020 adaptive Frangi, their MATLAB p-code not used): css_common
  .radius_profile = in-plane half-width at the 2-D centre line (EDT - half a voxel; 1/3/5-voxel bands -> 0.25/0.75/1.25
  mm at 0.5 mm). Columns cand_radius_med/p95_mm, cand_radius_cv; net_radius_med/p95_mm on the dark + vessel voxels
  of the connected network. PH3: lesions median 0.73 mm vs others 0.40 mm (phantom only). Quantised at the voxel size.
- Considered, not adopted now: LupoLab adaptive Frangi (protected .p files, erodes 17 px of brain surface, MRA-
  validated); brain-MRI-SSL encoders (2-D ResNet-50, no SWI in pretraining, no licence file) -> optional
  embed_candidates.py for the ranking model once ~10 expert-labelled patients exist.
- v4.18 from the external-repo review (BrainIAC perturbation test, ABIDE FreeSurfer morphometry): robustness.py
  re-runs detection on perturbed copies (ID_R<k>, removed afterwards): copy, in-plane / slice resampling (anti-
  aliased), bias field, noise, gamma, 1-voxel image shift against the labels; reports candidates, found again (centroid
  <= 5 mm, one-to-one), new, Spearman of matched scores, top-k overlap, known lesions kept -> review/ID_robustness.csv.
  PH3 (phantom): gamma 1 identical to the copy; known 6/6 kept under every perturbation; a 1-voxel shift raised
  candidates 7 -> 25 and slice 2 mm 7 -> 12 (registration / protocol sensitivity to check on P006).
  surface_features adds surf_sulc, surf_curv, surf_thickness_mm, surf_sulcal_frac (?h.sulc/curv/thickness).
  Repos reviewed and not adopted: nnU-Net cSS detector (needs many outlined positives), BrainIAC model (96^3,
  research-only licence), tumour toolkits (100-voxel cleanup, active contours, slice-level splits).
- v4.18 public controls: user's "6 new patients" (OpenNeuro ds004265) were MACAQUES (FUS safety study) - not used.
  OpenNeuro ds000221 (LEMON, CC0): 112 healthy adults with acq-mag_GRE + acq-phase_GRE (Siemens Verio 3 T, TE 17.3 ms,
  0.8 mm iso, LAI = axial, S-I axis 2), MP2RAGE T1w 1 mm, 3-D FLAIR; 26 aged 55-75 with GRE + 3-D FLAIR chosen.
  make_swi.py (homodyne sigma 4 mm per axial slice, paramagnetic sign auto from dark-on-magnitude voxels, mask^4;
  tests/check_make_swi.py: both sign conventions found, veins x0.03-0.04, calcium x1.00). run_lemon.sh: download
  (S3 curl), make_swi -> run_css --t1 --flair --recon --phase -> study_meta --split dev, ids P201+ in
  labels/id_map.csv, review/lemon_summary.csv. Presumed negative, not expert-read. MP2RAGE + recon-all untested.
- v4.18 from the microbleed papers (MixMicrobleedNet, MixMicrobleed, CenSynCMB 2026 - read in full): css_common
  .susceptibility_signal (dipole kernel in k-space, B0 = scanner z through the affine, phase 2 pi gamma B0 TE field,
  linear-gradient intravoxel dephasing |sinc|); make_synthetic --physics (cSS +depth ppm, veins +0.45, calcium -1.0,
  inside signal exp(-TE dR2*), SWI mask^4 on paramagnetic phase; --calcium K, --b0, --te-ms; work/HS_dphase);
  stress_test --physics. check_physics: cylinders match analytic field (parallel chi/3, perpendicular -chi/6), sphere
  +2:-1, B0 from affine; PH1 blooming TE10 4661 -> TE40 19390 voxels; calcium phase negative; 8/8 found (intensity
  model 7/8). evaluate.py: patient-grouped evaluation from labels/ (consensus or one reader; shown + excluded objects
  + deduped missed outlines; FP excludes uncertain; top-k; FROC; patient bootstrap CIs; --by strata; default split
  test); tests/check_evaluate.py with hand-computed answers. AIBL SWI (CSIRO 50304, 370 older adults, 3 T Siemens,
  0.93x0.93x1.75 mm, TE 20 ms) named as possible SWI control source - title says synthetic microbleeds: verify first.
- LEMON P201 on the user's Mac: make_swi sign -1 (median -0.28 rad, 72 %, n=17771; 0.4 % voxels darkened 42 %);
  recon-all on the MP2RAGE UNI aborted in topology fixing ('extra large defect', dura in wm.mgz: noisy background).
  v4.18 fix: mp2rage_robust.py (O'Brien 2014 robust combination from magnitude UNI+INV1+INV2, INV1 sign from UNI,
  lambda auto 1-32 until the border is < 10 % of the head; check_mp2rage: background SD 1182 -> 255, brain r 1.000);
  run_lemon downloads inv-1/inv-2, uses the robust T1, falls back to SynthSeg-on-T1 + mri_coreg if recon-all still
  fails; prep_anat moves an unfinished subjects/ID aside (ID_failed_<date>) before recon-all -i.
- v4.18.2: on P201 the border-based noise estimate read 10 (LEMON images are zero-padded/defaced at the border) ->
  lambda 1, speckle left around skull and neck (user's freeview). Noise now from INV2 in the band 5-20 mm outside an
  Otsu head mask; lambda 1..128 until < 5 % of that band is brighter than 25 % of the brain, stopping before the brain
  correlation drops below 0.95. check_mp2rage now zero-pads the border: old version lambda 1, speckle 83 -> 81 %; new
  lambda 32, 83 -> 0 %, brain r 1.000, WM median -6 %.
- v4.18.3 (LEMON P201 finished with recon-all, 6 candidates / 7 excluded, but QC 'pial edge agreement not
  measurable', edge_n 0, registration cost null):
  (1) css_qc boundary sample was |pial_sd| <= 0.75 x largest voxel; the voxel layers either side of a label boundary
  sit at >= the smallest voxel size, so on isotropic voxels (0.8 mm: 0.6 < 0.8) it was EMPTY. Now
  max(0.75 x largest, 1.01 x smallest) - unchanged for anisotropic scans. check_qc test 6 (PH3 relabelled 0.8 mm iso):
  old edge_n 0 / not measurable, new edge_n 4000, shift 0.02 mm, ok.
  (2) FreeSurfer 8 bbregister writes ID_swi2t1.dat.mincost (not .lta.mincost) -> cost was never read; both names now.
  (3) P203 stopped silently after '[0/4] importing + skull-stripping' on both routes (set -e, SynthStrip output sent
  to /dev/null): run_css now reports the failing tool and keeps work/ID_synthstrip.log. Cause still to be read there.

- v4.19 dimensionality on the surface (surface_features.py, recorded only): in one slice a vein and a cSS band are
  both lines; on the pial surface a vein is a curve and cSS a patch. New columns surf_width_mm (core projected area /
  geodesic length), surf_cross_mm (2 x p95 distance from the diameter path), surf_inradius_mm (largest distance from
  the footprint edge), surf_aspect (cross / length). run_css.sh now runs surface_features.py after detection when
  recon-all surfaces exist (failure only warns). PH6: cSS width 4.98-7.61 mm vs surface vein / sinus 0.47-0.93 mm;
  inradius >= 5.93 vs <= 1.93. Resolution limit (bands cut from lesion 5, 1/2/4/8 slices of 1.2 mm): width 0.65 /
  1.34 / 2.88 / 3.09 mm - a band only 1-2 slices high measures like a vein, so the feature can only help for cSS
  covering >= ~4 slices (~5 mm). The phantom cSS are wide by construction; real cSS widths are unknown until labelled.

- v4.19.1: the 2-at-a-time LEMON batch (8 threads each) crashed the user's laptop; only P201 finished. run_lemon now
  runs 1 subject at a time, CSS_THREADS (default 4, --threads N) for recon-all / SynthSeg, under nice -n 15.

- LEMON P201 with v4.19.1 (user's Mac): QC fix works on real data (label shift 0.22 mm, bias -0.09, scatter 1.51,
  edge measured); bbregister cost 0.96 read (> 0.8 -> low_confidence although the edge fit says aligned; limit not
  changed from one scan). 6 candidates in this presumed-healthy scan; user's visual read: 1, 2, 4, 5, 6 dark cortex /
  artifact, none veins (pre/postcentral, cuneus, precuneus). Surface width 1.6-5.6 mm, subpial: patch-shaped like cSS,
  so surf_width does NOT separate them. freeview default window shows this SWI white: use :grayscale=p1,p99.
- v4.19.2 regional + contralateral reference (detect_css.py, recorded only): region_z (candidate median vs its own
  cortical label's median, label MAD units), region_dark_frac (share of the label's cortex outside the candidate at
  least half as dark as the candidate, from the whole-cortex median), mirror_region_z (label vs the same label in the
  other hemisphere), cand_mirror_z (candidate vs that mirror label). tests/check_region.py: PH3 + a gyrus darkened
  25 % with its own label: cSS region_z -6.9..-9.2, dark frac <= 0.02; bilateral dark gyrus -0.45 / 0.78 / mirror 0;
  left-only dark gyrus -0.09 / 0.46 / mirror -3.6. Limits: real physiological darkness covers only part of a label
  (posterior precentral bank), receive-coil bias can make hemispheres differ, and cSS can be bilateral.
  Real data (v4.19.2, user's Mac): P006 known cSS (6/6 kept) region_z -3.64, -4.68, -4.37, -3.21, but the two in
  POSTCENTRAL cortex -2.35 / -2.18 (dark frac 0.25 each) - the same range as P201's normal motor-strip candidates
  (-1.43 / -1.79, 0.31 / 0.41) and P201 1/6 (-2.6). P201 3/5 at -3.65/-3.61 sit inside the cSS range. A -3 cut would
  lose 2 of 6 known lesions: region_z must NOT become a rule. Likely reason: cSS in naturally dark cortex, and
  extensive cSS darkens its own label. mirror_region_z ~0 everywhere (-1.0..+1.0), also in P006: uninformative so far.

- LEMON P202 / P204 (v4.19.1 batch, 1 at a time, 4 threads): both recon-all OK; 15 / 16 candidates (P201: 6); label
  shift 0.11 / 0.12 mm; QC low_confidence only from bbregister cost 0.964 / 0.954 (P201 0.962). False positives by
  region, 3 scans: inferior parietal 6, cuneus/pericalcarine 4, superior parietal 2, temporal 4, pre/postcentral 2 ...
  region_z medians -2.64 / -2.74 / -3.18 (33-56 % below -3) - inside the P006 known-cSS range, as expected.
- v4.19.3 css_qc: bbregister cost > 0.8 lowers the status only when the edge fit could not be measured; otherwise it
  is a note (qc["notes"]). check_qc test 7: cost 0.96 + measured edge fit -> ok with the note.

- LEMON GRE header (checked on P202): sform_code 2, qform 0; world box x -167..-1, y -195..9, z -162..-34 vs the
  MP2RAGE's x -85..90, y -131..108, z -122..133 (same session) - the GRE is NOT in scanner space, so freeview shows
  nothing when GRE and T1 are overlaid. Axis DIRECTIONS are correct: in SWI world coordinates of the registered labels,
  cerebellum below cerebrum, frontal in front of occipital, left on the left (P201/P202/P204 all True). Detection only
  uses relative positions; registration (mri_coreg init) copes. The right inferior parietal false positive sits at the
  same voxel position in P201/P202/P204 (44-46, 64-65, 43-48); user's read: dark cortex. Cause not established.

## Known issues / decisions
- SWI-only route (SynthSeg on SWI) gives a coarse pial surface → pial features less reliable.
- v4 ranking weights are hand-set from the definition, not trained; review all candidates.
- Medial cSS is mildly down-ranked by the midline factor (0.7).
- QSM/SEPIA not used (site decision). Phase only as a reading aid (calcium mimic); it cannot separate
  hemosiderin from deoxy-Hb.
- minIP series must not be used as input; the reviewer computes its own minIP slab for reading.
- Workbook (tools/build_workbook.py) builds the Excel validation workbook; recalculated with LibreOffice.

## Next steps (operating procedure: RUNBOOK.md)
1. P006 with --recon: check_known 6/6, QC file, surface_features (agreement < 2 mm); P011 count (expect 8).
2. Expert labels with fine codes (VS!) + missed-lesion outlines + study_meta split for every new case.
3. Patient-grouped evaluation script reading labels/ (FP per negative scan, top-k sensitivity, FROC, CIs).
4. Regularised ranking model with nested patient-level CV; review tiers from it; hard rules only after test.
5. Decide the phenotype scope for inferior temporal/occipital cSS (excluded by E1/E3/E8 by design).

## How to test
`bash tests/run_tests.sh` — v3 phantom (radial sulci, fissure, veins) end to end, fails if synthetic
sensitivity < 4/8; v4 phantom (check_v4.py: sensitivity, score_v4 AUC ≥0.8, ICH Dice, sulcal 3/4).
tests/check_stage2.py + check_consensus.py (v4.17), check_surfvein.py + check_overlay.py (v4.16, PH6), check_labels.py (v4.15) and check_qc.py (all called by run_tests.sh): v4.14 QC status, column contract, dropped log, extent at the
volume edge, label-shift monotonicity, transposed-orientation detection.
`bash tests/test_shell.sh` — run_css.sh / prep_anat.sh / run_all.sh flow with stubbed FreeSurfer.
No FreeSurfer or patient data needed.

## Conventions
- Python 3.11, numpy/scipy/nibabel only (+matplotlib for the reviewer); no network calls.
- Every script: CLI with `sys.argv`, reads `$CSS_BASE`, prints a one-line summary first.
- Prefer small, verifiable changes; run tests/run_tests.sh after any change to detect/score logic.
