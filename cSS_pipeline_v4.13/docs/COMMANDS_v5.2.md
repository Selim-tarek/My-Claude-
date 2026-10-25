# cSS pipeline v5.2 — all commands, in order, with recommendations

Covers every version up to v5.2 (v4 detector, v4.14–v4.19 QC / labels / surface / evaluation, v5.0 mimic layer and
ranker, v5.1.x fates / pieces / features / tiers, v5.2 patches and classifier). Details: RUNBOOK.md.
**Research software, not clinically validated. A human reads every candidate.** Study codes only (P012), never
names or MRNs; never commit patient data.

```zsh
setopt interactivecomments               # FIRST: otherwise zsh treats every '# ...' below as arguments
echo 'setopt interactivecomments' >> ~/.zshrc   # once, makes it permanent
S=~/css_project/scripts                  # set once per Terminal
```
Sections 0 and 1 can be pasted as they are. **Section 2 is a template**: set `P` to a real study code and replace
every `<...>` with a real path before running it, one step at a time - never paste it as a block.

## 0. Install / update (once per new version)

```zsh
cd ~/Downloads && unzip -o cSS_pipeline_v5.2.zip && cd cSS_pipeline_v5.2
conda env update -f environment.yml      # first time: conda env create -f environment.yml
zsh install.sh                           # backs up the old scripts, copies to ~/css_project/scripts
# open a NEW Terminal
conda activate css
bash tests/run_tests.sh | tail -3        # MUST end with: ALL TESTS PASSED
bash tests/test_shell.sh | tail -1       # MUST end: SHELL TESTS PASSED
```

## 1. Before you use any new switch: real-data regression (once per version)

```zsh
OLD=~/css_project/review/v3/P006_candidates.nii.gz; IDS=3,4,5,13,18,19
python $S/detect_css.py P006 -2.5 85 3 | head -4
python $S/check_known.py P006 --old $OLD --ids $IDS                          # MUST: 6 of 6 shown
CSS_USE_MIMIC=1    python $S/detect_css.py P006 -2.5 85 3 | head -4; python $S/check_known.py P006 --old $OLD --ids $IDS
CSS_MIXED_REVIEW=1 python $S/detect_css.py P006 -2.5 85 3 | head -4; python $S/check_known.py P006 --old $OLD --ids $IDS
CSS_USE_MIMIC=1    python $S/detect_css.py P201 -2.5 85 3 | grep "mimic layer"   # healthy: <= 1-2 with mimic_score > 0.5
python $S/detect_css.py P006 -2.5 85 3 > /dev/null; python $S/detect_css.py P201 -2.5 85 3 > /dev/null   # back to default
```
A switch becomes routine only if its line keeps 6/6. Otherwise leave it off and send the printed lines.

## 2. Per patient (template - edit, then run step by step)

```zsh
P=P0XX                                        # the real study code
# 2.1 split FIRST, before anyone looks at the output (dev = may tune on it; test = frozen pipeline only)
python $S/study_meta.py $P --split test

# 2.2 DICOM -> NIfTI (never the minIP)
R=~/css_project/raw/$P
dcm2niix -z y -f "%s_%d" -o $R/SWI   "<SWI DICOM folder>"
dcm2niix -z y -f "%s_%d" -o $R/T1    "<3-D T1 DICOM folder>"
dcm2niix -z y -f "%s_%d" -o $R/FLAIR "<FLAIR DICOM folder>"
grep -h '"ImageType"' $R/SWI/*.json          # pick the SWI (or magnitude) series, not minIP / phase

# 2.3 detection (recommended full route; recon-all 1-3 h)
export CSS_THREADS=4                          # one patient at a time on the laptop
CSS_QC_STRICT=1 CSS_MIXED_REVIEW=1 run_css.sh $P $R/SWI/<swi>.nii.gz \
    --t1 $R/T1/<t1>.nii.gz --flair $R/FLAIR/<flair>.nii.gz --recon \
    --phase $R/SWI/<phase>.nii.gz             # --phase only if the phase series exists
#   batch, no viewer: NOVIEW=1 run_all.sh P012 P013 P014   (logs/)
#   exit 3 = QC insufficient_quality under CSS_QC_STRICT: fix the input, do not read

# 2.4 QC - do not skip
python -c "import json; q=json.load(open('$HOME/css_project/review/${P}_qc.json')); print(q['status'], q.get('reasons'), q.get('notes'))"
freeview -v ~/css_project/data/${P}_swi.nii:grayscale=p1,p99 ~/css_project/work/${P}_t1seg_swispace.nii.gz:colormap=lut:opacity=0.3
#   check CORONAL and SAGITTAL too (an up-down slip is invisible axially). Off by > 1-2 mm:
#   prep_anat.sh $P --t1 ~/css_project/data/${P}_t1.nii.gz --recon && NOVIEW=1 run_css.sh $P
python $S/audit_case.py $P | tail -20        # inputs fresh, every object accounted for

# 2.5 (optional) lobar ICH: draw work/${P}_ich.nii.gz in freeview, then
python $S/detect_css.py $P

# 2.6 surface (recon-all cases; surface_features runs automatically in run_css.sh)
python $S/surface_overlay.py $P              # stops if surface-to-SWI fit > 2 mm -> back to 2.4

# 2.7 expert reading = the reference standard (blinded; excluded objects mixed in)
python $S/expert_sheet.py $P --include-excluded       # send review/${P}_expert_sheet.pdf, KEEP the _key.csv
python $S/expert_import.py $P expertA --letters "1C:5 2VS:4 3VX 4A 5U ..." --require-confidence
python $S/expert_import.py $P expertB --letters "..." --require-confidence
python $S/consensus.py --readers expertA,expertB --subjects $P
#   fill labels/adjudication_todo.csv, then: consensus.py --readers expertA,expertB --apply labels/adjudication_todo.csv

# 2.8 missed cSS + patient-level reference (needed for sensitivity)
python $S/label_missed.py $P expertA                  # from labels/${P}_missed_expertA.nii.gz
python $S/study_meta.py $P --status pos --ref-L 1 --ref-R 2 --motion 0 --reader expertA

# 2.9 score (from the expert calls; or your own reading with review_css.py)
python $S/expert_import.py $P expertA --letters "..." --score
python $S/review_css.py $P                   # default order (score_v4)
python $S/review_css.py $P --by-priority     # or: tier order (v5.1.3)
python $S/score_css.py $P --explain          # shows which sulci were linked and why

# 2.10 v5.2 patches (research; needs recon-all or T1 labels)
python $S/surface_aligned_patch.py $P
python $S/classify_surface_patches.py features $P
```

## 3. Across patients

```zsh
python $S/feature_report.py                                   # features vs calls
python $S/evaluate.py --split test --source consensus --by field_T   # sensitivity, FP/scan, FROC, CIs
python $S/agreement.py --a ~/css_project/review_session1 --b ~/css_project/review_session2   # intra/inter-rater
python $S/robustness.py P006 --old ~/css_project/review/v3/P006_candidates.nii.gz --ids 3,4,5,13,18,19
python $S/stress_test.py --hosts P011 --physics --calcium 2 --tag v52phys      # dev hosts only

# learned ranking - only when the label counts below are reached
python $S/train_ranker.py --split dev --source consensus      # -> review/ranker_model.json + ranker_model_card.md
python $S/rank_candidates.py P020                             # run_css.sh does this for dev studies automatically
python $S/classify_surface_patches.py train --split dev --source consensus      # ablations A-G
python $S/classify_surface_patches.py predict P020
```

Public healthy controls (LEMON): `run_lemon.sh ~/Documents/lemon --download --subjects "sub-010239"` then
`--summary` (one subject at a time).

## 4. Recommendations

| Setting | Recommendation | Why |
|---|---|---|
| T1 + `--recon` | **Always** | pial surface, sulcal scoring, surface steps; without T1 every case is low confidence |
| FLAIR | Yes when acquired | acute cSAH / thrombosis context (soft penalty only, never an exclusion) |
| `--phase` | Yes when the phase series exists | calcification mimic (E10); otherwise the phase rule is inert |
| `CSS_QC_STRICT=1` | **On** for study runs | a failed case stops (exit 3) instead of producing unreadable candidates |
| `CSS_MIXED_REVIEW=1` | On, once section 1 shows 6/6 | can only ADD objects back for review, never removes one |
| `CSS_USE_MIMIC=1` | Only after section 1: P006 6/6 and P201 quiet | E10 / E11 remove objects; v5.0.0 E11 once removed a true lesion (fixed in v5.1.0) |
| `review_css.py --by-priority` | Optional reading aid | tier cut-offs are fitted to phantoms; still read every candidate |
| `train_ranker.py` | Wait for >= 10 labelled dev patients (refuses < 3) | below that the intervals are too wide to beat score_v4 |
| `classify_surface_patches.py train` | Wait for >= 5 labelled dev studies with cSS AND veins | refuses fewer; start from the SWI-only ablation (A), not G_all |
| Any model as the default order | Only after it beats score_v4 on **test** studies (`evaluate.py`) | never tune on test |
| `CSS_THREADS=4`, one patient at a time | Yes on the laptop | 2 parallel runs at 8 threads crashed it (v4.19.1) |
| `CSS_KEEP_ALL`, `CSS_SATO_LEGACY`, `CSS_ALLOW_ORIENT` | Debugging only | disable exclusions / change filters / bypass the orientation stop |
| Expert codes | Use the specific vein codes (`VS`, `VC`, `VT`, `VX`) and a confidence `:1-5` | they are the hard negatives both models need |
| Phenotype scope | Decide before the study | inferior temporal / occipital cSS over the tentorium is excluded by design (E1 / E3 / E8) |

Re-run `check_known.py P006 ...` after every update; the 6 known lesions must stay shown.
