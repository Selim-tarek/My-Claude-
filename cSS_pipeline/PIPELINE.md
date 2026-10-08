# cSS pipeline v4 — step-by-step guide

Semi-automatic detection and grading of **cortical superficial siderosis (cSS)** on SWI, with
T1 (and optionally FLAIR) anatomy. Research software, not for clinical use.

---

## Published definitions this pipeline follows

- **cSS** (Charidimou et al., *Neurology* 2017;89:2128, after the consensus criteria): "well-defined,
  homogeneous hypointense curvilinear signal loss (black) on T2*-GRE or SWI **outlining the outer
  surface of cerebral cortex, within the adjacent subarachnoid space, or both**". It is limited to the
  convexity sulci and spares the brainstem, cerebellum and spinal cord (Charidimou et al., *Brain* 2015).
- **Multifocality scale** (same 2017 paper), per hemisphere: 0 none; 1 = one sulcus or up to 3
  immediately adjacent sulci; 2 = 2 or more non-adjacent sulci or more than 3 adjacent sulci. Total
  0–4. More than 3 adjacent sulci count as multifocal because that much involvement implies more
  than one bleeding source. Interrater weighted κ 0.87; rated blinded.
- **ICH rule** (same paper): cSS "contiguous or potentially anatomically connected with any lobar
  ICH" is not scored. Every cSS area must be separated from any lobar ICH by **at least 3 unaffected
  sulci**, or **at least 2 unaffected sulci (at multiple axial levels) if the haematoma has no
  superficial communication path along the convexity in 3-D**. v4.3 counts these sulci automatically
  (`score_css.py`, needs `--recon`):
  - fewer than 2 sulci in between → excluded, if you drew the ICH mask; only a warning with the
    automatic mask;
  - exactly 2 sulci in between → kept, but flagged for you to check the two extra conditions.
- **Minimum standards** (Charidimou, *AJNR* 2016;37:E43). Of the minimum criteria, these are used
  here:
  - **blooming** on T2*/SWI;
  - **differentiation from microbleeds**: microbleeds are small (generally 2–5 mm), well-defined,
    round/oval, and at least half surrounded by brain parenchyma. v4.4 flags candidates ≤10 mm with
    ≥50 % parenchyma around them as `cmb_like`. Beware "microbleed pearls", rows of tiny
    microbleeds along the surface;
  - **T1 or FLAIR for anatomical confirmation** of the gyral location;
  - **no infratentorial siderosis**;
  - **exclude mimics**: vessel flow voids, thrombosed vessels, petechial haemorrhagic
    transformation of infarcts, calcium deposits;
  - FLAIR/PD-bright subarachnoid space → call it **acute cSAH**;
  - **record location and number of sulci**, and **record other lesions within 1 cm**
    (microbleeds, DWI lesions).
- **SWI physics** (Vaccarino et al., *Bioengineering* 2025;12:473; Barbosa et al., *Radiol Bras*
  2015;48:93):
  - Magnitude alone cannot tell paramagnetic material (haemosiderin, deoxy-Hb in veins) from
    diamagnetic material (**calcium**). The **phase** can.
  - The phase sign convention depends on the vendor: Siemens/Canon are "left-handed", GE/Philips
    "right-handed". So always compare with a vein on the same image: calcium has the **opposite
    sign** to veins and blood products.
  - Echo time and the SWI phase-mask choice change the image, and calcium can look dark like
    haemosiderin.
  - v4.5:
    - `run_css.sh ... --phase <phase.nii>` adds a **filtered-phase column** to the reviewer;
    - the score file records **TE, field strength, vendor and voxel size**, because blooming
      (apparent cSS size) depends on them;
    - minIP is a venography tool, which is why it is a reading aid here and never the input.
- **Boston criteria v2.0** (Charidimou et al., *Lancet Neurol* 2022;21:714) count cSS **by gyri**
  for the CAA diagnosis:
  - a single focus, even if it extends to a second adjacent gyrus, = 1 haemorrhagic lesion;
  - multifocal cSS (gyri separated by uninvolved areas, or ≥3 adjacent gyri) = ≥2 lesions.
  - Probable CAA = ≥2 strictly lobar haemorrhagic lesions (ICH, microbleeds, cSAH/cSS), or 1 lobar
    haemorrhagic lesion plus a white-matter feature (>20 CSO perivascular spaces in one hemisphere,
    or multispot WMH), with a compatible presentation, age ≥50 and no other cause.
  - v4.4 reports `boston2_css_lesions` (0 / 1 / 2 = ≥2). This is only the cSS part, not a diagnosis.
  - Interrater κ for multifocal cSS was 0.94.
- **SWI vs T2\*-GRE:** SWI gives higher multifocality scores than GRE, so the sequence is now saved
  with every score (`sequence` column) and the two must not be pooled.
- **Quantification** (van Harten et al., *NeuroImage: Clinical* 2023;38:103447, Leiden):
  - Method: 2-D Frangi vesselness (σ = 1 voxel), rater-placed seeds, 3-D 6-neighbour region
    growing, and false positives **removed by a human rater**.
  - They report that the filter catches all cSS "but also several other hypointense tubular
    structures such as some veins". The vein problem is the known open problem, and a human
    review step is the published standard.
  - Volume reproducibility: ICC 0.995. Voxel overlap between sessions: Dice 0.75.
  - They propose **% of sulci affected** as less sequence-dependent than volume. v4.3 reports it
    per hemisphere (`L_sulci_pct`, `R_sulci_pct`).

## 0. What changed in v4, and why

v3 asked: *"is this dark, line-shaped and near cortex?"* Veins answer "yes" too. That is why, on
P006, the v3 ranking score had an AUC of only 0.53 against your reader calls.

v4 tests each part of the **textbook definition** separately (STRIVE-2, Charidimou multifocality
scale, Radiopaedia):

| Definition of cSS | What v4 measures (CSV column) | cSS | Vein / mimic |
|---|---|---|---|
| Hypointensity **outlining the cortical surface** (subpial / subarachnoid) | `pial_dist_mm`: signed distance to the pial surface (− = inside cortex, + = in CSF); `bank_frac` | ≈ −1 to 0 mm, bank_frac high | middle of the sulcal CSF (> +0.5 mm) |
| **Curvilinear**, following gyral contours | `tube_ratio`: 3-D Hessian shape, 0 = flat sheet, 1 = tube | < 0.4 (sheet) | > 0.5 (tube) |
| Follows the gyrus | `surface_alignment`: the sheet lies parallel to the cortex | ≈ 1 | lower |
| **Tram-track**: both banks of a sulcus | `tram_frac`: walks across the sulcus and tests the opposite bank, with brighter CSF between them | > 0 | ≈ 0 (one central line) |
| Not part of the venous tree | `vein_tree_mm`: extent of the connected dark-tube network (**recorded only**: on real SWI the network merges brain-wide) | short | long, branching |
| Asymmetric | `mirror_dark_frac`: darkness at the mirror point in the other hemisphere | low | high (normal veins are paired) |
| **Remote from lobar ICH** | `ich_dist_mm`, `near_ich_suggest` (≤ 5 mm) | far | near → reader presses **i** |
| **Supratentorial** | `infratentorial` | 0 | 1 → classical superficial siderosis, reported separately |
| Chronic: **no FLAIR hyperintensity** | `flair_csf_z`, `flair_bright` (only with FLAIR) | not bright | bright → acute convexity SAH |
| Mimic: **cortical vein thrombosis** | `flair_ctx_z`, `flair_ctx_bright`: FLAIR signal of the cortex next to it | not bright | bright cortex → check for thrombosis |
| Counted **by sulci** (focal 1–3, disseminated > 3) | `score_css.py` with Destrieux labels | | |

Also changed in v4:
- **Search zone = the pial band** (cortex + subarachnoid CSF). White matter and deep nuclei are
  excluded. v3 dilated the cortex into white matter, which is darker than grey matter on SWI.
- **Local darkness reference excludes dark voxels** (two passes), so extensive cSS no longer
  darkens its own reference.
- **Elongation filter 2.0** (was 2.5). A tram-track lesion is two parallel plates, elongation
  about 2.4, and v3 dropped those. Round microbleeds stay about 1–1.5.
- **T1 anatomy** (`prep_anat.sh`): accurate pial surface, plus Destrieux sulci when recon-all is
  available.
- **Sulcal scoring** (`score_css.py`): "adjacent sulci" means sulci that touch each other or
  border the same gyrus. It replaces the 10 mm Euclidean rule that turned P006 left from 1 into 2.
- **Bug fixes**:
  - the batch script now runs detection;
  - reader calls are tied to the lesion itself, not its rank;
  - every candidate is reviewed by default, and unreviewed candidates are reported;
  - hemisphere assignment is fixed;
  - a stale SynthSeg segmentation is no longer reused;
  - the synthetic test requires ≥30 % lesion coverage.
- **Reviewer** shows an 8 mm minIP slab and the definition features for every candidate.

> **Status: honest.** On both phantoms all definition features separate cSS from veins.
> **The v4 ranking weights are NOT validated on real patients yet.** They are transparent
> rules, not a trained model. That is why every candidate is now reviewed, so the score does
> not depend on the ranking. Validate on P006 and the next cases before trusting the ranking
> (step 6).

---

## 1. One-time setup (Mac)

```zsh
cd ~/Downloads/cSS_pipeline          # this folder
zsh install.sh                       # copies scripts to ~/css_project/scripts (backs up old ones)
# open a NEW Terminal window, then:
conda activate css
conda install -n css -c conda-forge matplotlib scikit-image -y     # if not already there
bash tests/run_tests.sh              # phantom tests, ~3 min, must end with ALL TESTS PASSED
bash tests/test_shell.sh             # checks run_css.sh / prep_anat.sh flow (FreeSurfer stubbed)
```

FreeSurfer 8.2 must be set up in your shell (`$FREESURFER_HOME`, licence file). `recon-all`
is only needed for sulcal scoring.

---

## 2. Which series to export

| Series | Use? | Why |
|---|---|---|
| **SWI** (processed SWI image, or the SWI magnitude) | **YES, main input** | Detection runs on this |
| **3-D T1** (MPRAGE / SPGR / BRAVO, ~1 mm) | **YES**, strongly recommended | Accurate pial surface, sulci, hemisphere |
| **FLAIR** (3-D or 2-D) | Yes, if available | Separates chronic cSS from acute convexity SAH |
| **SWI minIP** | **NO, not as input** | Slab projection: it merges veins and cSS from neighbouring slices and destroys the 3-D shape analysis. The reviewer computes its own minIP panel from the SWI. |
| SWI phase / filtered phase | Optional, reading aid only (`--phase`) | Both haemosiderin and deoxy-Hb are paramagnetic, so phase does not separate cSS from veins, but it **does** flag calcium (opposite sign), an official cSS mimic |

Convert DICOM, one folder per series:

```zsh
P=P007; R=~/css_project/raw/$P
dcm2niix -z y -f "%s_%d" -o $R/SWI   "<SWI DICOM folder>"
dcm2niix -z y -f "%s_%d" -o $R/T1    "<T1 DICOM folder>"
dcm2niix -z y -f "%s_%d" -o $R/FLAIR "<FLAIR DICOM folder>"
ls $R/*                              # file names = <series number>_<series description>
grep -h '"ImageType"' $R/SWI/*.json  # skip any file whose ImageType contains MIN_IP / MNIP / PHASE / "P"
```

Use **study codes only** (P007), never names or MRNs. Never commit or share patient data.

---

## 3. Run one patient

### Route A — fast (minutes): SWI + T1 + FLAIR, SynthSeg on the T1
```zsh
run_css.sh P007 ~/css_project/raw/P007/SWI/<swi>.nii.gz \
           --t1 ~/css_project/raw/P007/T1/<t1>.nii.gz \
           --flair ~/css_project/raw/P007/FLAIR/<flair>.nii.gz
```
Steps that run:
1. SWI: import, then SynthStrip (CSF kept, so cSS at the surface is not stripped).
2. `prep_anat.sh`:
   - T1: import, SynthStrip, `mri_synthseg --parc --robust`.
   - Register SWI → T1 with `mri_coreg` (rigid).
   - Labels → SWI grid with `mri_vol2vol --inv --interp nearest`.
   - FLAIR: import, SynthStrip, `mri_coreg` → SWI, resample.
3. `align_seg.py`: uses the T1 labels and masks them to the SWI coverage.
4. `detect_css.py`: candidates and features → `review/P007_candidates.csv/.nii.gz`.
5. freeview opens with the candidates (unless `NOVIEW=1`).

Scoring in Route A uses the **Euclidean approximation** (no sulcal labels).

### Route B — best (1–3 h): adds recon-all → sulcal scoring
```zsh
run_css.sh P007 <swi> --t1 <t1> --flair <flair> --recon
```
Same as Route A, except:
- `recon-all -all` runs on the T1 (log: `work/P007_recon.log`);
- registration uses `bbregister --t2 --init-coreg` (boundary-based, more accurate);
- Destrieux labels (`aparc.a2009s+aseg`) are brought onto the SWI grid, so `score_css.py`
  scores **by sulci**.

Tip: start recon-all for the next patients overnight, e.g.
`recon-all -s P008 -i <t1> -all -threads 8 &`. `prep_anat.sh` re-uses a finished recon-all
automatically.

### Route C — SWI only (no T1 available)
```zsh
run_css.sh P007 <swi>
```
Uses SynthSeg on the SWI (coarser cortex, so the pial-distance features are less reliable).

### Batch (no viewer)
```zsh
run_all.sh P007 P008 P009        # each needs data/<ID>_swi.nii; logs in ~/css_project/logs/
```

---

## 4. Quality control (2 minutes, do not skip)

`prep_anat.sh` prints a freeview command such as:
```zsh
freeview -v data/P007_swi.nii work/P007_t1seg_swispace.nii.gz:colormap=lut:opacity=0.3 \
         work/P007_flair_swispace.nii.gz:visible=0 work/P007_a2009s_swispace.nii.gz:colormap=lut:opacity=0.3:visible=0
```
Check in three planes:
- the coloured cortex ribbon follows the SWI cortex, especially at the vertex and in the
  occipital lobes;
- the FLAIR edges line up when you toggle it on and off.

**If it is misaligned** by more than about 1–2 mm, delete `work/P007_t1seg_swispace.nii.gz` and
re-run `run_css.sh P007`. It falls back to SWI SynthSeg. Note the case in the workbook.

**Lobar ICH present?** The detector finds large, compact dark haematomas automatically
(`work/P007_ich_used.nii.gz`; check it). To override, draw the haematoma (including its
hemosiderin rim) in freeview as a new volume, save it as `work/P007_ich.nii.gz`, and re-run
`python ~/css_project/scripts/detect_css.py P007`.

---

## 5. Review and score

```zsh
python ~/css_project/scripts/review_css.py P007          # ALL candidates (default)
```
Screen layout:
- **Top left:** the whole slice.
- **Bottom left:** an **8 mm minIP slab**. Veins become continuous, branching tubes; cSS stays a
  band along the cortex.
- **Right:** three adjacent slices, raw on top and outlined underneath.
- **Title:** region, volume, darkness, and the **definition features** (pial distance,
  plate/tube shape, follows-surface, tram-track %, vein tree, mirror darkness), plus flags:
  vein-like, ICH distance, infratentorial, FLAIR-bright.

Keys: **y** cSS · **v** vein · **o** normal · **a** artifact · **i** near ICH · **u** unsure ·
**b** back · **[ ]** slice · **q** finish. Every key press is saved, so you can quit and resume.

Reading checklist for each candidate (radiological criteria):
1. Is it **on the cortical surface** (subpial or in the subarachnoid space), following the gyral
   contour over consecutive slices? → cSS-compatible.
2. **Tram-track**, i.e. both banks of the same sulcus dark? → strongly favours cSS.
3. In the minIP panel, does it **continue as a tube** into a branching network or towards the
   sinus? Is it a **dot** on perpendicular slices? → vein.
4. Is the same structure present at the **mirror location**? → favours a normal vein, or normal
   iron in the motor cortex.
5. **Contiguous with a lobar ICH**? → press **i** (reported, not scored).
6. **Infratentorial** (cerebellum, brainstem)? → not cSS. Think of classical superficial
   siderosis.
7. **FLAIR-bright sulcus** at the same place? → acute convexity SAH, not chronic cSS.

When you press **q**, scoring runs automatically (`mark_css.py` → `score_css.py`):
```
P007  cSS multifocality score: 3/4   (disseminated)   [sulcal (Destrieux)]
  left : score 1  foci 2  clusters 1  sulci 2
  right: score 2  foci 3  clusters 3  sulci 4
  sulci: ctx_lh_S_postcentral, ctx_lh_S_intrapariet_and_P_trans, ...
  volume: candidates ... mm3, grown (full extent) ... mm3
```
- **Per hemisphere:** 0 = none; 1 = one sulcus, or ≤ 3 immediately adjacent sulci; 2 = ≥ 2
  non-adjacent sulci, or > 3 sulci. **Total 0–4.**
- **STRIVE-2:** focal = 1–3 sulci; disseminated = > 3 sulci.
- `WARNING: N candidates were not reviewed`: the score is a **lower bound**. Finish the review.
- Results are written to `review/P007_score.json` and `review/css_scores.csv` (one row per
  patient), and the grown cSS mask to `review/P007_css_mask.nii.gz`.

Enter the result in the workbook (Case Log). **Review before you look at the expert score**
(blinding).

---

## 5a. Fewer false positives (v4.8)

Reader feedback on P006 (duplicates, specks, obvious veins, skull-base artifact, faint normal cortex)
led to these steps. **Excluded candidates are never deleted.** They are listed with the reason in
`review/ID_excluded.csv` (and `ID_excluded.nii.gz` to view in freeview).

| Problem | Fix |
|---|---|
| One dark line reported 2–3 times | Pieces closer than **3 mm** are merged into one candidate (this also joins the two banks of a tram-track), **but only while the merged object stays ≤ 40 mm long**; longer chains (on P006: veins chained into 10–17 cm objects) keep their pieces separate |
| Tiny specks | Candidates smaller than **6 mm** overall are excluded (curvilinear cSS is longer) |
| Microbleeds | Small and ≥ half surrounded by brain tissue (AJNR 2016 rule) → excluded |
| Obvious veins | **Tubular** (tube_ratio > 0.45), **out in the middle of the sulcal CSF** (> 0.5 mm outside the pial surface) and **no tram-track** → excluded |
| Skull-base artifact | Artifact-zone labels (orbitofrontal, temporal pole, …) in the lowest 30 % of the brain → excluded |
| Faint normal dark cortex | Less than **0.75 SD** darker than its 1–2 mm surroundings ("well-defined" in the definition) → excluded |

| Cerebellum / tentorium lines (v4.9) | ≥ 30 % of the candidate within 5 mm of cerebellum or brainstem → excluded (cSS is supratentorial by definition; infratentorial siderosis is a different disease) |
| Basal cisterns (v4.9) | ≥ 30 % within 5 mm of brainstem, ventral diencephalon, hippocampus or amygdala → excluded (circle of Willis, posterior cerebral artery, basal vein of Rosenthal) |
| Skull base (v4.9) | Lowest 20 % of the **cerebrum's** height (any label), or artifact-zone labels in the lowest 40 % → excluded |

**v4.9.2: when the position rules apply.** They backfired on P010 (a partial slab, no T1): synthetic
lesions were excluded as "skull base" or "basal cisterns".
- **Skull-base and basal-cistern rules** run only with **T1 anatomy** and when the scan actually
  **includes the skull base**. A slab cut off at the bottom has no skull base in it.
- **Without a T1**, the cerebellum rule needs ≥ 50 % of the candidate near the cerebellum.
- **Microbleed-like** is now a flag in the reviewer, not an exclusion.
- The **merge gap** grows with slice thickness (2 × slice, 3–4 mm). The **elongation filter**
  applies only to single-piece candidates.
- The first output line states the anatomy source and the coverage.

- Safety checks: on both phantoms every cSS lesion is still kept (6/6 and 8/8), while half the
  phantom veins are removed.
- On a real scan, check that your **known** lesions survive:
  ```
  python ~/css_project/scripts/check_known.py P006 --old ~/css_project/review/v3/P006_candidates.nii.gz --ids 3,4,5,13,18,19
  ```
- To switch every exclusion off (e.g. to compare): `CSS_KEEP_ALL=1 python ~/css_project/scripts/detect_css.py P006`

## 5b. Expert reading (the reference standard)

Your own calls are a "non-expert + tool" reading. The **reference** must come from an expert.
Make it quick for them:

1. **Make the blinded picture sheet:**
   ```
   python ~/css_project/scripts/expert_sheet.py P006
   ```
   This creates `~/css_project/review/P006_expert_sheet.pdf`:
   - one row per candidate (SWI zoom | outlined | 8 mm minIP);
   - **random order** with new numbers, and **no scores**, so the expert isn't biased;
   - an instruction page and an answer grid.

   The key (sheet number → candidate) is saved separately in `P006_expert_key.csv`. **Don't send
   the key**, only the PDF. It contains no names, only the study code.
2. The expert writes a letter for every number: **C** cSS · **V** vein · **N** normal/other ·
   **H** contiguous with an ICH · **U** unsure.
3. **Import the answers.** Numbers you don't list count as N:
   ```
   python ~/css_project/scripts/expert_import.py P006 expertA --letters "1V 2V 3C 4N 5U ..."
   ```
   or, shorter:
   ```
   python ~/css_project/scripts/expert_import.py P006 expertA --css 3,17,40 --unsure 5
   ```
   Add `--score` to compute the 0–4 score from the expert's calls. Your own calls are not changed.
4. For a second reader, use the same sheet with another name (`expertB`). This gives inter-rater
   agreement.
5. `feature_report.py` now uses the expert calls automatically when they exist. Use
   `--source reader` to use yours instead.

## 5c. Synthetic lesions (v4.10)

`make_synthetic.py` inserts artificial cSS into a cSS-negative scan, so we know exactly where the
lesions are. The v3 generator makes them **realistic**: the P006 stress test showed that the old
lesions were too big and blob-like, so size and branching separated them "for the wrong reason".

- **Thin curvilinear lines** traced along the cortex–CSF boundary: 1 voxel wide, sometimes 2
  (blooming), **5–30 mm** long, on 2–5 slices.
- **Tram-track** in 40 % of lesions: the opposite bank of the same sulcus gets a matching line.
- **Vein decoys** (default 4 per synthetic scan in the stress test): dark tubes 0.5–0.9 mm in
  radius and 15–40 mm long, running through the sulcal CSF. The stress test reports each feature's
  AUC **vs these decoy veins** as well.
- `--legacy` reproduces the old generator for comparison.

```
python ~/css_project/scripts/stress_test.py --hosts P006 --depths 0.6,0.45,0.3 --seeds 1,2 --tag p006_v3
```

## 6. Validation plan (what makes this publishable)

1. **Re-run P006 with v4.** First back up the v3 review: `mkdir -p ~/css_project/review/v3 && cp ~/css_project/review/P006_* ~/css_project/review/v3/`.
   Then run `run_css.sh P006 <swi> --t1 <t1> --recon` and review again, ideally after a few days
   or blinded.
   - Target: left = 1 (the two adjacent foci are now scored by sulci), total 3.
   - Old v3 calls are ignored automatically after re-detection (the cand_ids change), so they
     can never be attached to the wrong lesion.
2. After **every** case: `python ~/css_project/scripts/feature_report.py`. It shows, for each
   feature, the AUC for cSS vs everything else across all reviewed cases.
   - Expectation from the definition: `tube_ratio`, `pial_dist_mm`, `vein_tree_mm` and
     `mirror_dark_frac` below 0.5; `tram_frac`, `bank_frac` and `surface_alignment` above 0.5.
   - Compare `score_v4` with `score_v3`. Treat a feature as useful only if its AUC is ≥ 0.75 or
     ≤ 0.25 consistently across patients.
3. With about **15 reviewed cases**, replace the hand-set `score_v4` with logistic regression on
   these features. Use **leave-one-patient-out** cross-validation, not pooled candidates.
4. **Inter-rater reliability:** a second reader scores a subset blinded; report κ for presence
   and weighted κ for the 0–4 score. Compare the pipeline + reader result with the expert
   reference (the workbook Validation sheet does this).
5. Include **cSS-negative CAA or age-matched scans** to measure the number of false-positive
   candidates per scan.
6. Synthetic stress test on real hosts (as before): `python ~/css_project/scripts/stress_test.py --tag v4`.

---

## 7. Method details (for the methods section)

1. **Pre-processing.**
   - SWI skull-stripped with SynthStrip, extra-axial CSF kept.
   - T1 skull-stripped and segmented, either by SynthSeg 2.0 (`--parc --robust`) or by
     recon-all (aparc+aseg, aparc.a2009s+aseg).
   - Rigid registration SWI → T1: `mri_coreg`, or `bbregister --t2` with recon-all.
   - Labels resampled to the SWI grid with nearest-neighbour interpolation.
   - FLAIR registered rigidly to the SWI.
2. **Pial surface.**
   - Tissue = labelled non-CSF voxels, holes filled (ventricles).
   - Signed Euclidean distance to the tissue boundary in mm.
   - Surface normal = gradient of the distance map smoothed at σ = 1 mm.
3. **Intensity normalisation.**
   - Clip at the 99.5th percentile.
   - z = (I − local cortical mean) / cortical MAD×1.4826.
   - Local mean is a Gaussian-weighted (σ = 10 mm) mean over cortex. Second pass excludes
     cortex with z < −2.
4. **Candidate generation.**
   - Search zone: −3.5 to +4 mm from the pial surface, minus white matter, deep grey nuclei and
     a 3 mm skull-strip rim.
   - Per-slice Sato ridge filter (σ 1–3 voxels, dark ridges).
   - Strong voxels: z < −2.5 and ridge > 85th percentile. Weak voxels: z < −2.0 and ridge > 80th
     percentile.
   - Hysteresis connected components per hemisphere.
   - Filters: ≥ 25 mm³, ≥ 2 slices, elongation ≥ 2.0.
5. **Features.**
   - Multi-scale 3-D Hessian (σ = 0.8 and 1.6 mm, scale-normalised); eigenvalues sorted by
     magnitude. `tube_ratio` = |λ2|/|λ3| at voxels with λ3 > 0. `surface_alignment` = |e3·n|.
   - Tram-track: rays from bank voxels along n up to 10 mm. The opposite bank is where the ray
     re-enters tissue. A tram is counted when the opposite bank has z < −2 and the sulcal midline
     is > 1 SD brighter than both banks.
   - Venous network: dark (z < −2) tubular (tube_ratio > 0.5) voxels at pial distance > −1.5 mm;
     component bounding-box diagonal.
   - Mirror: plane fitted (SVD) to the interhemispheric boundary, then reflection.
   - ICH: lobar tissue with z(WM) < −3, ≥ 500 mm³, inscribed radius ≥ 2.5 mm, median > 3 mm
     under the pial surface, not midline (< 5 mm), not at the skull-strip rim, not in artifact zones.
   - FLAIR: median robust z of CSF within 3 mm.
6. **Ranking score (v4, untrained).**
   - √darkness × (0.25 + on-surface) × (0.25 + sheetness) × (0.5 + alignment) × (1 + tram).
   - Multiplied by 0.5 for an artifact zone or a long structure, 0.7 for midline, 0.8 for
     mirror > 0.6, and 0.3 for infratentorial. (`vein_tree_mm` is recorded only: on P006 the
     dark-tube network merged brain-wide, 255 mm for most candidates.)
   - Darkness is square-rooted because veins were darker than cSS on P006.
7. **Scoring.**
   - Each accepted focus is assigned to Destrieux sulcal labels within 4 mm (≥ 15 % of voxels).
   - Two sulci are adjacent if they touch or border a common gyral label.
   - Charidimou 0–2 per hemisphere; STRIVE-2 focal/disseminated.
   - Without Destrieux: 3 mm foci and 10 mm clusters (approximation).
   - Volume: seeded region growing into dark voxels (z < −1.5) within 5 mm.

---

## 8. Troubleshooting

| Problem | Fix |
|---|---|
| `run_all.sh` shows FAILED | Read `~/css_project/logs/<ID>_run.log` |
| Candidates hug the wrong place / labels shifted | QC step 4, then delete `work/<ID>_t1seg_swispace.nii.gz` |
| "old-format calls IGNORED" | Expected after a v4 re-run; v3 calls pointed at v3 ranks. Review again. |
| Many candidates in orbitofrontal / temporal pole | Susceptibility artifact zones; they are flagged `artifact_zone` |
| Need the old ranking for comparison | `CSS_RANK=v3 python ~/css_project/scripts/detect_css.py P007` (both scores are always in the CSV) |
| recon-all failed | `work/<ID>_recon.log`; Route A still works |
