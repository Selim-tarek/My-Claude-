# cSS pipeline v5.0 — runbook

Step-by-step operation of the whole pipeline, from DICOM to a validated research score. Method details and
the reasoning behind each rule are in [PIPELINE.md](PIPELINE.md); the audit and plan in
[docs/AUDIT_v4.13_and_FP_plan.md](docs/AUDIT_v4.13_and_FP_plan.md).

> **Research software, not clinically validated.** The detector proposes candidates; a human reads every one.
> All ranking weights and exclusion thresholds are hand-set from phantoms and two development patients
> (P006, P011). None has been tested on held-out patients yet.

Replace `P012` with your study code. Never use names or MRNs. Never commit patient data.

---

## The pipeline at a glance

```
 0  install + self-tests                                   (once)
 1  assign the validation split                            study_meta.py --split        (before anything else)
 2  export DICOM -> NIfTI                                  dcm2niix
 3  run detection                                          run_css.sh  (strip, anatomy, detect)
 4  quality control  ── stop here if insufficient_quality  review/P012_qc.json + freeview
 5  (optional) lobar ICH mask
 6  (recon-all cases) surface features + surface view      surface_features.py, surface_overlay.py
 7  expert reading (the reference standard)                expert_sheet.py -> expert_import.py
 8  missed cSS + patient-level reference                   label_missed.py, study_meta.py
 9  score                                                  expert_import.py --score  (or review_css.py)
10  per-case bookkeeping                                   feature_report.py, surface_overlay --color-by expert
11  across cases: reliability, stress tests, evaluation    agreement.py, stress_test.py   (section 11)
13  v5.0: mimic layer, model ranking, strict QC             CSS_USE_MIMIC, train_ranker.py, CSS_QC_STRICT  (section 13)
```

---

## 0. One-time setup (Mac)

```zsh
cd ~/Downloads/cSS_pipeline_v4.16
conda env create -f environment.yml   # once (or: conda env update -f environment.yml)
zsh install.sh                     # copies scripts to ~/css_project/scripts, backs up old ones
# open a NEW Terminal window
conda activate css
bash tests/run_tests.sh            # ~10 min, must end with: ALL TESTS PASSED
bash tests/test_shell.sh           # must end with: SHELL TESTS PASSED
```
FreeSurfer 8.x must be set up (`$FREESURFER_HOME`, licence). `recon-all` is needed for sulcal scoring and for the
surface steps (6, 10). All scripts use `$CSS_BASE` (default `~/css_project`).

Re-run both test scripts after every update of the scripts.

---

## 1. Assign the validation split — before anyone looks at the output

```zsh
python ~/css_project/scripts/study_meta.py P012 --split test      # or: --split dev
```
- `test` cases are only ever scored by a frozen pipeline; nothing is tuned on them.
- The split cannot be changed afterwards without `--force` (and the change is logged).
- P006, P010, P011 and every case already used to design rules are `dev`.

---

## 2. Export the series

```zsh
P=P012; R=~/css_project/raw/$P
dcm2niix -z y -f "%s_%d" -o $R/SWI   "<SWI DICOM folder>"
dcm2niix -z y -f "%s_%d" -o $R/T1    "<3-D T1 DICOM folder>"
dcm2niix -z y -f "%s_%d" -o $R/FLAIR "<FLAIR DICOM folder>"
grep -h '"ImageType"' $R/SWI/*.json
```

| Series | Use |
|---|---|
| SWI magnitude **or** processed SWI | Main input. Pick one kind per study protocol and keep it (they have different contrast) |
| 3-D T1 | Strongly recommended (pial surface, sulci). Without it, the case is `low_confidence_review` |
| FLAIR | Recommended (acute convexity SAH, cortical oedema flags) |
| SWI phase | Optional, reading aid for calcium (`--phase`) |
| **minIP** | **Never** as input |

---

## 3. Run detection

```zsh
run_css.sh P012 ~/css_project/raw/P012/SWI/<swi>.nii.gz \
    --t1 ~/css_project/raw/P012/T1/<t1>.nii.gz \
    --flair ~/css_project/raw/P012/FLAIR/<flair>.nii.gz \
    --recon                                           # recon-all, 1-3 h: sulcal scoring + surface steps
#   [--phase ~/css_project/raw/P012/SWI/<phase>.nii.gz]
```
Several patients without the viewer: `run_all.sh P012 P013 P014` (logs in `~/css_project/logs/`).

What runs: SynthStrip (CSF kept) → `prep_anat.sh` (recon-all + bbregister, or SynthSeg on the T1 + mri_coreg;
FLAIR registration) → `align_seg.py` → `detect_css.py`.

The detector prints three lines:
```
P012: anatomy T1, coverage includes skull base, merge gap 3.0 mm
P012: 2908 raw components -> 41 candidates [23 excluded: 12 vein, ...] ...
P012: QC ok, pial edge agreement 0.6 mm, 512 components below the generation gates [review/P012_qc.json; reported only]
```

Main outputs:

| File | Content |
|---|---|
| `review/P012_candidates.csv/.nii.gz` | candidates to read (cand_id = rank) and all features |
| `review/P012_excluded.csv/.nii.gz` | objects removed by a rule, with `reason`, `rule_version` (E1–E9), `exclusion_confidence` |
| `work/P012_dropped.csv/.nii.gz` | components removed before feature extraction (size, slices, elongation, no cortex) |
| `review/P012_qc.json` | case QC (step 4) |
| `work/P012_vessels.nii.gz`, `work/P012_dark.nii.gz` | vessel map; dark-voxel map for region growing |
| `work/P012_reg.json` | registration method and cost |
| `review/P012_manifest.json` | provenance: input file hashes, parameters, package and FreeSurfer versions |
| `work/P012_gen.nii.gz` | per voxel, which detection stage it passed (used to explain missed cSS) |

---

## 4. Quality control (do not skip)

1. **Read the status** in the third detector line or `review/P012_qc.json`:
   - `insufficient_quality` → **stop**. The SWI is not stored with the slice axis last, or the labels cover
     too little of the brain. Fix the input; do not read the candidates.
   - `low_confidence_review` → continue, but anatomy-based features and exclusions are less reliable. Note the
     reasons (SWI-only anatomy, labels shifted > 1.5 mm against the SWI, label boundary > 1.5 mm inside/outside the
     SWI cortex edge on average, registration cost, slices > 2 mm).
   - The QC file separates `label_shift_mm` (rigid misregistration, with its direction `label_shift_vec_mm`),
     `edge_bias_mm` (labels systematically too small > 0 or too large < 0) and `edge_scatter_mm` (noise: veins,
     partial volume — not a registration problem).
   - `ok` → continue.
   The thresholds are provisional; nothing in the pipeline acts on the status yet.
2. **Look at the registration** with the command `prep_anat.sh` printed:
   ```zsh
   freeview -v data/P012_swi.nii work/P012_t1seg_swispace.nii.gz:colormap=lut:opacity=0.3
   ```
   The cortex ribbon must hug the SWI cortex in all three planes, especially at the vertex, occipital poles and
   temporal lobes. `edge_agreement_by_region_mm` in the QC file says where the fit is worst.
   **Always look at a coronal or sagittal view**: an up-down slip (third number of `label_shift_vec_mm`) is
   invisible in axial slices (P011: 3.9 mm up-down, axial view looked perfect).
   If it is off by more than ~1-2 mm, redo the registration first: re-run
   `prep_anat.sh P012 --t1 data/P012_t1.nii.gz --recon` (bbregister on the cortex boundary instead of
   `mri_coreg`; the imported T1 copy can be given as `--t1`), then `NOVIEW=1 run_css.sh P012`. Only if that also
   fails: delete `work/P012_t1seg_swispace.nii.gz` and re-run `run_css.sh P012` (falls back to SynthSeg on the
   SWI; the case becomes `low_confidence_review`).
   `prep_anat.sh` stops with an ERROR at the first failed step (e.g. a wrong `--t1` path); it never leaves the
   old labels in place silently.

---

## 5. (Optional) lobar ICH mask

The detector finds large lobar haematomas automatically (`work/P012_ich_used.nii.gz`; check it). To set it
yourself, draw the haematoma including its haemosiderin rim in freeview, save as `work/P012_ich.nii.gz`, then
`python ~/css_project/scripts/detect_css.py P012`. A drawn mask makes the ICH sulcal rule exclude
(rather than only warn about) cSS within 2 sulci of the ICH.

---

## 6. Surface steps (recon-all cases only)

```zsh
python ~/css_project/scripts/surface_features.py P012     # pial-surface features (run_css.sh does this automatically since v4.19 when recon-all surfaces exist); surf_width_mm / surf_inradius_mm: curve (vein) vs patch (cSS)
python ~/css_project/scripts/surface_overlay.py P012      # candidates (red) / excluded (blue) on the surface
```
- Both scripts build the surface→SWI transform from `subjects/P012/mri/orig.mgz` and `work/P012_swi2t1.lta`,
  print the median distance of the pial vertices to the SWI label boundary for both directions, and **stop if
  the better one is > 2 mm**. A stop means: go back to step 4.
- `surface_overlay.py` prints a `freeview -f ...inflated:annot=...:overlay=...` command. Click a vertex to read
  the object id (+ candidate, − excluded). `review/P012_surf_list_objects.csv` lists the Destrieux sulcus and
  gyrus under each object.
- Outputs are **recorded only**; they do not change candidates, ranks or scores.

---

## 7. Expert reading (the reference standard)

Your own calls (step 9, `review_css.py`) are a "non-expert + tool" reading. The reference comes from an expert,
blinded to scores and ranks.

1. **Make the blinded sheet**, with the excluded objects mixed in (unmarked):
   ```zsh
   python ~/css_project/scripts/expert_sheet.py P012 --include-excluded
   ```
   → `review/P012_expert_sheet.pdf` (send this) and `review/P012_expert_key.csv` (**do not send**).
   Each row: axial SWI, outlined, 8 mm minIP, coronal, sagittal.
2. **The expert writes one code per number**, optional confidence `:1`–`:5`:

   | Code | Meaning | Code | Meaning |
   |---|---|---|---|
   | `C` | cSS | `A` / `AM` | air-bone artifact / motion artifact |
   | `V` | vein (unspecified) | `N` / `NI` | normal or other / normal iron-rich cortex |
   | `VS` | **vein along the cortical surface** | `MB` / `CA` | microbleed / calcification |
   | `VC` | mid-sulcal vein | `SAH` | acute convexity SAH |
   | `VT` | transcortical vein | `LN` / `HI` | laminar necrosis / haemorrhagic infarct |
   | `VX` | vein at a sinus | `IS` | infratentorial siderosis |
   | `TV` | thrombosed vein | `H` / `U` | contiguous with lobar ICH / uncertain |

   Please use the specific vein codes (especially `VS`): they are the hard-negative set the next steps need.
3. **Import**:
   ```zsh
   python ~/css_project/scripts/expert_import.py P012 expertA --letters "1C:5 2VS:4 3VX 4A 5U ..."
   ```
   Unlisted numbers count as `N` (flagged `defaulted`). If the expert calls an **excluded** object `C`, the
   import warns: that exclusion rule removed real cSS — note it.
4. **Second reader** (same sheet): `expert_import.py P012 expertB --letters "..."`.
   **Blinded re-read** weeks later: add `--session 2`.
5. **Consensus**: `consensus.py --readers expertA,expertB --subjects P012` → objects both readers put in the
   same class get a consensus label; disagreements go to `labels/adjudication_todo.csv`. A third reader (or a
   consensus meeting) fills `adjudicated_label` (a code such as `C` or `VS`) and `adjudicator`, then
   `consensus.py --readers expertA,expertB --apply labels/adjudication_todo.csv`. Kappa (cSS vs not, six
   classes) with patient-bootstrap CIs → `labels/agreement_stats.json`. **The consensus is the reference for
   evaluation**, not either reader alone.

All answers go to `labels/candidate_labels.csv`.

---

## 8. Missed cSS and the patient-level reference

1. **Missed cSS.** The expert outlines every cSS focus that has **no number** on the sheet, in freeview on the SWI:
   `freeview -v data/P012_swi.nii review/P012_candidates.nii.gz:colormap=lut` → File › New Volume (template: the
   SWI) → draw → save as `labels/P012_missed_expertA.nii.gz`. Then:
   ```zsh
   python ~/css_project/scripts/label_missed.py P012 expertA
   ```
   For each outline it records whether it was excluded (which rule), dropped (which gate) or never generated,
   and if never generated, why: outside the search zone (skull-strip rim / white-matter label / beyond the pial
   band), not dark enough, no ridge response, or no strong seed → `labels/missed_lesions.csv`. **Without this step sensitivity cannot be measured.** If outlining is not
   practical, ask for slice and location and draw it yourself from those notes.
2. **Patient-level reference**:
   ```zsh
   python ~/css_project/scripts/study_meta.py P012 --status pos --ref-L 1 --ref-R 2 --motion 0 --reader expertA
   ```
   (`--status pos|neg|uncertain`; `--ref-L/R` = expert multifocality 0–2 per hemisphere.) Protocol and QC fields
   are filled from `review/P012_qc.json` → `labels/study_meta.csv`.

---

## 9. Score

From the expert's calls (the reference):
```zsh
python ~/css_project/scripts/expert_import.py P012 expertA --letters "..." --score
```
Or your own interactive reading of every candidate (keys `y` cSS, `v` vein, `o` normal, `a` artifact,
`i` near ICH, `u` unsure, `b` back, `[ ]` slice, `q` finish → scores automatically):
```zsh
python ~/css_project/scripts/review_css.py P012
```

Output (`review/P012_score.json`, one row per patient in `review/css_scores.csv`):
- cSS multifocality 0–4 (Charidimou 2017; per hemisphere 0 / 1 / 2), sulcal with Destrieux labels, else a
  Euclidean approximation; STRIVE-2 focal / disseminated;
- Boston v2.0 cSS lesion count (by gyri), ICH sulcal rule, infratentorial siderosis reported separately;
- foci, sulci, volume, grown volume, cortical surface area covered and % per hemisphere;
- sequence, field, TE, voxel size (sizes are comparable only within one protocol).

A `WARNING: N candidates were not reviewed` line means the score is a lower bound.

`--by-tier` (v4.18, optional) shows tier-1 candidates first. `review_tier` = 2 when the candidate's own flags
point to a mimic (artifact zone, vein-like, microbleed-like, vein evidence > cSS evidence; reason in
`tier_reason`). It is NOT validated: every candidate is still shown and must still be read, and expert sheets
keep their own order. Use it only after `check_known.py` shows the known lesions in tier 1.

---

## 10. Per-case bookkeeping

```zsh
python ~/css_project/scripts/feature_report.py                       # features vs calls, all reviewed cases
python ~/css_project/scripts/surface_overlay.py P012 --color-by expert   # expert calls on the surface
```
After **any** detector change, re-check the known lesions of earlier cases (v4.18 also prints each kept
lesion's review tier and what its two ends do: stops / runs_on / branches / midline):
```zsh
python ~/css_project/scripts/check_known.py P006 --old ~/css_project/review/v3/P006_candidates.nii.gz --ids 3,4,5,13,18,19
```

---

## 11. Across cases

| Task | Command |
|---|---|
| Reader agreement on candidate calls | `consensus.py --readers expertA,expertB` (step 7.5) |
| Inter-/intra-rater agreement of two scoring sessions (incl. weighted kappa of the 0–4 score) | `cp -R review review_session1` (later `review_session2`); `agreement.py --a review_session1 --b review_session2 --col surface_mm2 --plot results/agreement.png` |
| Volume / surface vs multifocality (ceiling effect) | `agreement.py --vs-score ~/css_project/review --col surface_mm2 --plot results/vs_score.png` |
| Robustness of one scan to protocol / image changes (v4.18; ~8 detection runs) | `robustness.py P006 --old review/v3/P006_candidates.nii.gz --ids 3,4,5,13,18,19` |
| **Evaluation against expert labels (v4.18)**: patient-grouped sensitivity (all / top-k), FP per scan and per negative scan, FROC, bootstrap CIs, strata | `evaluate.py --split test --by field_T,slice_group` (dev split: `--split dev`, optimistic) |
| Synthetic stress test with susceptibility physics (v4.18: blooming, calcium mimics) | `stress_test.py --hosts P011 --physics --calcium 2 --tag v418phys` |
| Synthetic stress test on cSS-negative hosts (dev only) | `stress_test.py --hosts P011 --depths 0.6,0.45,0.3 --seeds 1,2 --surface-veins 2 --tag v416` |
| Shortlist exclusion rules (dev only, hypothesis generation) | `rule_test.py --pos P006 --old ... --ids ... --neg P011` |

**Rules for using the results**
- Thresholds, weights and rules are chosen on `dev` cases only, with all candidates of a patient kept together.
- Report false positives per **negative** scan, sensitivity within the top 5/10/20 candidates, and per-protocol
  results, with patient-bootstrap confidence intervals — not only AUC.
- A new hard exclusion is adopted only if, on `test` cases, it removes **no** adjudicated cSS and reduces false
  positives per negative scan.

**Phenotype scope (a study decision, not a setting):** candidates on the inferior temporal / occipital surface
over the tentorium and near the skull base are excluded by design (rules E1, E3, E8; all listed in
`ID_excluded.csv`). Decide before the study whether such cSS is outside the phenotype, or must be retained
(then those rules have to become review flags), and state it in the methods.

**Not built yet (planned, in this order):** a patient-grouped evaluation script that reads `labels/`; a
regularised ranking model validated by patient-level cross-validation; review tiers from that model. Until
then, every candidate is read and the hand-set `score_v4` only orders the list.

---

## 12. Public control scans: OpenNeuro ds000221 (MPI-Leipzig LEMON), v4.18

Healthy volunteers (CC0), Siemens 3 T, gradient-echo magnitude + phase at 0.8 mm isotropic (TE 17.3 ms), MP2RAGE T1,
3-D FLAIR. Use: **presumed cSS-negative controls, split dev** - false positives per scan, robustness to a new
scanner/protocol. Never sensitivity, never test. `make_swi.py` builds the SWI from magnitude + phase (homodyne
high-pass, paramagnetic phase mask^4, phase sign read from the image and printed - check once that veins get darker).

```zsh
zsh ~/css_project/scripts/run_lemon.sh ~/Documents/lemon --download --subjects "sub-010239" --parallel 1
```
First ONE subject: check `logs/P201_run.log`, `review/P201_qc.json` (label_shift_mm, status) and freeview in
coronal view (MP2RAGE background noise can disturb recon-all skull stripping). Then all 26 (overnight):
```zsh
zsh ~/css_project/scripts/run_lemon.sh ~/Documents/lemon --download --parallel 2
zsh ~/css_project/scripts/run_lemon.sh ~/Documents/lemon --summary
```
IDs P201... are fixed in `labels/id_map.csv`; every subject is entered in `labels/study_meta.csv` as split dev with
the note "PRESUMED cSS-negative, not expert-read". Have an expert read a sample before calling any candidate a false
positive.

---

## 13. v5.0: mimic layer, trained ranking, strict QC (all opt-in)

With no switch set, v5.0 gives the same candidates, exclusions and order as v4.19; the new columns are added at
the end of every candidates / excluded CSV.

| Switch / script | What changes | Default |
|---|---|---|
| `CSS_USE_MIMIC=1` | score x0.5 with FLAIR evidence (acute cSAH / vein thrombosis); exclusions **E10** calcification (phase), **E11** microbleed cluster | off: mimic columns recorded only |
| `CSS_QC_STRICT=1` | QC `insufficient_quality` -> no candidates, "case cannot be scored" (exit 3); `low_confidence_review` -> reviewer banner, model score x0.8; stale T1 labels stop `align_seg.py` | off |
| `train_ranker.py` / `rank_candidates.py` | `score_model`, `rank_model` columns; `review_css.py --by-model` | no model until you train one |

The phase rule needs the SWI **phase** imported with the scan (`run_css.sh ... --phase <phase.nii.gz>` ->
`data/P012_phase.nii.gz`); without it the rule is inert and the FLAIR and microbleed rules still run.

### 13a. Before any switch becomes routine: real-data regression (run on the Mac)

```zsh
S=~/css_project/scripts
# 1. defaults: the known lesions are still shown
python $S/detect_css.py P006 -2.5 85 3 | head -4
python $S/check_known.py P006 --old ~/css_project/review/v3/P006_candidates.nii.gz --ids 3,4,5,13,18,19
#    MUST: 6 of 6 known lesions still shown for review
# 2. the new columns are there, and quiet on the known lesions
python -c "import pandas as pd; d=pd.read_csv('$HOME/css_project/review/P006_candidates.csv').set_index('cand_id'); print(d.loc[[6,10,18,19,33,44],['mimic_score','mimic_class','phase_z','mb_cluster_n']])"
#    (use the candidate numbers check_known printed in step 1)
# 3. mimic layer on: still 6 of 6
CSS_USE_MIMIC=1 python $S/detect_css.py P006 -2.5 85 3 | head -4
python $S/check_known.py P006 --old ~/css_project/review/v3/P006_candidates.nii.gz --ids 3,4,5,13,18,19
# 4. a healthy scan: at most 1-2 candidates with mimic_score > 0.5 (none is expected)
CSS_USE_MIMIC=1 python $S/detect_css.py P201 -2.5 85 3 | grep "mimic layer"
python -c "import pandas as pd; d=pd.read_csv('$HOME/css_project/review/P201_candidates.csv'); print(d[d.mimic_score>0.5][['cand_id','region','mimic_class','mimic_reason']])"
# afterwards: back to the default output
python $S/detect_css.py P006 -2.5 85 3 > /dev/null; python $S/detect_css.py P201 -2.5 85 3 > /dev/null
```
If any step fails, keep the switch off and send the printed lines: they name the candidate and the rule.

### 13a2. Per-case audit (v5.1)

`python $S/audit_case.py P012` prints, from the existing outputs: input hashes and geometry (orientation, obliquity,
S-I axis), acquisition, whether T1 labels / FLAIR / phase / Destrieux / surface features are present and fresh,
registration and QC, and what happened to every generated object (not generated, dropped per gate, excluded per rule,
shown) - the counts must add up (`review/P012_audit.json`).

### 13b. Expert labels -> trained ranking -> apply

1. Expert reading as in section 7 (sheets, `expert_import.py`, `consensus.py`), on **dev-split** studies.
   At least 3 labelled patients; with fewer than ~10 the numbers are noisy.
2. `python $S/train_ranker.py --split dev --source consensus` -> `review/ranker_model.json`. It prints the
   leave-one-patient-out AUC of the model **and of score_v4 on the same objects**, with patient-bootstrap
   intervals. Use the model only if it beats score_v4 clearly (intervals, not point estimates).
3. `python $S/rank_candidates.py P006 --allow-training-study` and compare `rank_model` with `cand_id` for the
   known lesions. From then on `run_css.sh` ranks **dev** studies automatically; `review_css.py P012 --by-model`
   shows that order (`--by-legacy` / default = the score_v4 order).
4. Test-split studies are never trained on (`train_ranker.py` refuses `--split test`) and are not ranked
   automatically: rank them by hand with the frozen model, then `evaluate.py --split test`.

### 13c. Mixed merged objects (v5.1.1, opt-in)

`CSS_MIXED_REVIEW=1 python $S/detect_css.py P012 -2.5 85 3` keeps an object whose merged pieces are part vein-like,
part cSS-like for review instead of letting a vein rule (E5 / E6 / E7) exclude it (`mixed_object_review` = 1). Without
the switch, `mixed_object` is recorded only. Check on P006 that 6/6 known lesions stay shown with the switch on.

### 13d. Review-priority tiers (v5.1.3)

`python $S/review_css.py P012 --by-priority` reads tier 1 (no strong mimic evidence) first, then 2 (uncertain), 3
(likely vein, convergent evidence), 4 (microbleed-like / artifact zone); tier 5 = case QC insufficient. The tier is in
`review_priority_rank` / `review_priority_reason`. It is a reading order only - read every candidate. Check on P006:
the 6 known lesions should be tier 1 or 2; if any is tier 3, send its reason.

### 13e. Surface-aligned patches and the patch classifier (v5.2, research)

```zsh
python $S/detect_css.py P012 -2.5 85 3                 # writes work/P012_pial_sd.nii.gz
python $S/surface_aligned_patch.py P012                 # review/P012_patches/*.npz + P012_patches.csv
python $S/classify_surface_patches.py features P012     # review/P012_patch_features.csv
# once >= 5 dev studies carry expert labels (cSS vs veins):
python $S/classify_surface_patches.py train --split dev --source consensus      # ablations A-G, model JSON
python $S/classify_surface_patches.py predict P020      # never a training study
python $S/train_ranker.py --split dev --source consensus   # now nested LOPO + review/ranker_model_card.md
```
Read the ablation table: a channel is worth keeping only if it raises the LOPO AUC beyond the interval of the
SWI-only model. Neither model changes the reading order until it beats score_v4 on test studies (evaluate.py).

---

## Where things are (`$CSS_BASE`)

```
raw/P012/{SWI,T1,FLAIR}/      dcm2niix output
data/P012_swi.nii             skull-stripped SWI (input to everything)
subjects/P012/                recon-all
work/                         labels on the SWI grid, registrations, vessel/dark maps, dropped log, reg.json
review/                       candidates, excluded, QC, surface outputs, expert sheet + key, calls, scores
labels/                       candidate_labels.csv, missed_lesions.csv, study_meta.csv, missed-lesion masks
results/                      stress tests, agreement plots
logs/                         run_all.sh logs
```

## What remains unvalidated

Every exclusion threshold and the ranking (hand-set; E7 chosen in-sample on P006/P011); the QC status
thresholds; the vein-tracking and surface features (phantom regression only); sensitivity to faint cSS and to
cSS over the skull base/tentorium (excluded by design); behaviour on 1.5 T, GRE, other vendors and thick
slices; the FreeSurfer LTA handling on real recon-all output (checked automatically against the labels at every
run). The multifocality score itself is a published, validated visual scale; the pipeline's version of it has
not yet been compared with expert scores beyond one patient.
