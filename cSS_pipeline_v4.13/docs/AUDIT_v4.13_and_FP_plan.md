# cSS pipeline v4.13: code audit, evidence review and false-positive reduction plan

Status: **audit and plan only. No pipeline file has been changed.** Implementation waits for approval
(deliverable 7). This is research software. Nothing here is clinically validated, and nothing in this
document should be read as evidence that it is.

Contents
1. Code audit
2. Literature evidence table
3. Ranked false-positive reduction plan
4. Data schema for expert labels and hard negatives
5. Leakage-safe validation protocol
6. Staged implementation plan
7. Test-suite gaps and additions
8. What remains unvalidated

Baseline used for this audit: the unchanged v4.13 package, `bash tests/run_tests.sh` run in a clean
Python venv (numpy 2.5, scipy 1.18, scikit-image 0.26; no FreeSurfer). Results:
`bash tests/run_tests.sh` → **ALL TESTS PASSED**; `bash tests/test_shell.sh` → **SHELL TESTS PASSED**. The log
excerpts are quoted where they matter below.

---

## 1. Code audit

### 1.1 Candidate-generation path (voxel → candidate), `scripts/detect_css.py`

| Step | Lines | What happens |
|---|---|---|
| Inputs | 96-104 | `data/ID_swi.nii` (skull-stripped SWI; processed SWI *or* magnitude, see 1.6 #10), `work/ID_seg_swispace.nii.gz`. `T1_ANAT` = `work/ID_t1seg_swispace.nii.gz` exists and is newer than the SWI |
| Masks | 109-114 | `brain` = filled `I>0`; `cortex` = labels ≥1000 or 3/42; `edge` = 2 mm in-plane rim; `surface` = outside brain or label 24; `near_surface` = 2 mm in-plane dilation of `surface` |
| Pial signed distance | 116-123 | `tissue` = filled (brain ∧ seg>0 ∧ seg≠24). `pial_sd` = EDT(¬tissue) − EDT(tissue), in mm (anisotropic `sampling=vox`: correct). Normals = gradient of `pial_sd` smoothed at σ=1 mm, in mm units |
| Darkness z | 125-139 | Clip at the 99.5th percentile. Global cortical median/MAD → `sd`. Local mean = Gaussian(σ=10 mm) of cortex intensities. Pass 2 recomputes the local mean excluding cortex with z<−2. `z = (Ic − local)/sd` (local mean, **global** SD) |
| Search zone | 141-147 | `shell` = brain ∧ −3.5 ≤ pial_sd ≤ +4 mm ∧ ¬(WM ∪ deep), plus `edge`. `zone` = shell ∧ ¬rim (rim = 3-D EDT to the mask edge < 3 mm) |
| Ridge | 149-156 | 2-D Sato (σ = 1, 2, 3 **pixels**) on each array slice `[:, :, k]`. `strong` = zone ∧ z<−2.5 ∧ R>85th percentile of R in shell. `weak` = zone ∧ z<−2.0 ∧ R>80th percentile |
| Hysteresis | 158-174 | Hemisphere map from nearest L/R label. 26-connected components of `weak`, per hemisphere, kept if they contain a `strong` voxel |
| Merge | 176-203 | Components within the merge gap `min(max(3, 2·vox[2]), 4)` mm are merged per hemisphere. A merged group whose bounding-box diagonal is >40 mm is split back into its pieces |
| Size/shape gates | 382-399 | **Silently dropped** (not written anywhere): volume <25 mm³, <2 array slices, elongation <2.0 (single-piece only), and no cortical label within 2 voxels unless infratentorial |
| Features | 400-547 | See 1.3 |
| Exclusions | 561-604 | Rule cascade; excluded rows → `review/ID_excluded.csv/.nii.gz` with reason |
| Output | 605-631 | `review/ID_candidates.csv/.nii.gz` sorted by `score` (v4 by default). Also `work/ID_dark.nii.gz` (growing mask), `work/ID_vessels.nii.gz`, `work/ID_ich_used.nii.gz` |

### 1.2 Rules that remove or demote a candidate

| # | Rule | Threshold | Type | Needs T1 | Written to excluded.csv? |
|---|---|---|---|---|---|
| G1 | volume | <25 mm³ | hard, generation | no | **no (silent)** |
| G2 | slices | <2 array slices | hard, generation | no | **no (silent)** |
| G3 | elongation (single-piece only) | <2.0 | hard, generation | no | **no (silent)** |
| G4 | no cortex label within 2 voxels and not infratentorial | — | hard, generation | no | **no (silent)** |
| E1 | infratentorial / tentorial | `infratentorial` flag, or `infra_frac` ≥0.3 (T1) / ≥0.5 (SWI-only) within 5 mm of cerebellum/brainstem/4th ventricle | hard exclusion | no | yes |
| E2 | basal cisterns | `basal_frac` ≥0.3 within 5 mm of brainstem, ventral DC, hippocampus, amygdala | hard exclusion | yes | yes |
| E3 | skull base | `rel_height` <0.2 of cerebrum height, and coverage includes the skull base | hard exclusion | yes | yes |
| E4 | speck | `extent_mm` <6 | hard exclusion | no | yes |
| E5 | vein in sulcal CSF | `pial_dist_mm` >0.5 ∧ `tube_ratio` >0.45 ∧ ¬(`tram_frac` ≥0.2) | hard exclusion | no | yes |
| E6 | vessel into WM | `vessel_wm_mm` ≥4 ∧ `axis_normal` ≥0.5 ∧ ¬(`tram_frac` ≥0.4) | hard exclusion | no | yes |
| E7 | off the cortex | `cortex_frac` <0.30 ∧ `pial_dist_mm` >0.30 ∧ ¬(`tram_frac` ≥0.2) | hard exclusion | yes | yes |
| E8 | skull-base artifact zone | `artifact_zone` ∧ `rel_height` <0.4, coverage includes skull base | hard exclusion | no (!) | yes |
| E9 | faint / ill-defined | `edge_contrast` <0.75 | hard exclusion | no | yes |
| F1 | `cmb_like` | extent ≤10 mm ∧ parenchyma ≥0.5 | reviewer flag | — | — |
| F2 | `near_ich_suggest` | ICH distance ≤5 mm | reviewer flag | — | — |
| F3 | `flair_bright`, `flair_ctx_bright` | robust z >3 | reviewer flag | — | — |
| F4 | `vein_like` | pial >0.3 ∧ tube >0.4 | reviewer flag | — | — |
| R1 | score_v4 multipliers | artifact ×0.5, midline ×0.7, long (>40 mm in slice direction) ×0.5, infratentorial ×0.3 | ranking | — | — |

Notes:
- E8 runs without T1 anatomy, even though E2/E3 were restricted to T1 anatomy in v4.9.2 after P010. On the SWI-only route `artifact_zone` comes from SynthSeg-on-SWI labels, the route `CLAUDE.md` calls "too coarse".
- E1, E3 and E8 remove true cSS on the inferior temporal and basal occipital surface. `CLAUDE.md` v4.9 says this is a known trade-off. It changes the phenotype the tool can detect, so it should be stated as a scope limit in any methods text.
- The `CSS_KEEP_ALL=1` switch disables E1-E9 but not G1-G4.

### 1.3 Features (CSV columns), meaning and likely failure modes

| Column | Intended meaning | Likely failure modes |
|---|---|---|
| `volume_mm3`, `n_slices`, `elongation`, `extent_mm` | size and shape | Depend on voxel size, slice thickness and blooming (TE, field). `n_slices` counts array axis 2, whatever its anatomical direction. `extent_mm` is the bounding-box diagonal (not path length) and is **underestimated at the low edge of the volume** (1.6 #B1) |
| `darkness_z` | mean −z | z uses a global cortical MAD: depends on TE, field, processed SWI vs magnitude, receive-coil bias. Veins were darker than cSS on P006 (AUC 0.37) |
| `surface_contact` | fraction within 2 mm (in-plane) of CSF/outside | Depends on how sulcal CSF is labelled (1.6 #B2): route-dependent |
| `cortex_frac`, `cortex_dist_mm` | in/near cortex | Registration error (1-2 mm moves the whole value); NN label resampling onto thick SWI slices; atrophy |
| `pial_dist_mm`, `bank_frac` | signed distance to pial boundary | Voxel-level Euclidean distance to a *label* boundary, not to a pial surface. Quantised to the grid (the phantom gives −1.13 for every cSS). Thresholds of 0.3 and 0.5 mm are smaller than one voxel on most protocols (P011 was a 2.5 mm reformat) |
| `tube_ratio` | 3-D Hessian λ2/λ3: plate (0) vs tube (1) | Scale-limited to σ 0.8/1.6 mm: a 1-voxel cSS line on 2-mm slices looks tubular; thick slices make all structures plate-like in the slice direction. Computed on `Ic`, which is 0 outside the mask (1.6 #B3) |
| `surface_alignment` | \|e3 · n\|: sheet parallel to cortex | Also ≈1 for a **tube lying on the surface** (one of its two cross-section directions is the normal). Cannot separate surface-parallel veins from cSS by construction |
| `tram_frac` | opposite bank dark, brighter CSF between | Single normal ray per voxel. Fails for oblique/curved sulci, thick slices, closed sulci (no CSF between banks), and single-bank cSS (valid phenotype). Requires ≥5 hits else NaN |
| `vein_tree_mm` | extent of connected dark tubular network | Saturates brain-wide on real SWI (255 mm). Recorded only |
| `mirror_dark_frac` | darkness at reflected point | Venous anatomy is asymmetric; cSS can be bilateral. Reflection uses voxel coordinates and a fitted plane, so it is sensitive to head tilt. `z` is not masked outside the brain (1.6 #B4) |
| `ich_dist_mm`, `near_ich_suggest` | proximity to lobar ICH | Auto-ICH misfired on P006 (skull base, sinus) |
| `infratentorial`, `infra_frac`, `basal_frac`, `rel_height` | position | Label quality; `rel_height` assumes a full cerebrum in the FOV (gated by coverage check) |
| `artifact_zone`, `midline_zone` | DK region names / ≤5 mm from midplane | Region names come from the majority label, which may be a neighbour in thick slices |
| `branch_per10mm` | skeleton branch points | 2-D per array slice; length uses `vox[0]` for every step (diagonals undercounted); noisy for short objects |
| `surface_gradient` | outer vs inner cortex darkness | Needs ≥5 voxels each; meaningless at 2.5 mm voxels |
| `parenchyma_frac`, `cmb_like` | AJNR microbleed rule | Same `tissue` definition issue as `pial_dist` |
| `flair_csf_z`, `flair_bright`, `flair_ctx_z`, `flair_ctx_bright` | acute cSAH / cortical oedema | CSF reference = label 24 only (1.6 #B2); FLAIR rigidly registered; no inflow/flow-artifact handling |
| `edge_contrast` | "well-defined" | Uses the same z map; low SNR or blooming reduces it; P006 known cSS min 2.17, so E9 at 0.75 is loose |
| `vessel_ext_mm`, `vessel_wm_mm`, `vessel_run_mm`, `axis_normal` | vessel continuation, orientation | Bounded to 20 mm but still on a thresholded tube map; `vessel_run_mm` uses the candidate's straight long axis, so it fails for curved veins; `axis_normal` is high for cSS wrapping crown→fundus |
| `css_evidence`, `vein_evidence` | rule-based 0-1 summaries | Not probabilities, untrained |
| `score_v3`, `score_v4`, `score` | ranking | Hand-set; P006 AUC ~0.5 against non-expert calls |

### 1.4 Dependence on resolution, vendor, field, TE, orientation and registration

| Calculation | Resolution | Vendor/field/TE/processing | Orientation | T1↔SWI registration |
|---|---|---|---|---|
| Sato σ = 1,2,3 **pixels** | yes: 0.47-1.4 mm at 0.47 mm, 1.2-3.6 mm at 1.2 mm | — | 2-D in array slices | — |
| z thresholds −2.5/−2.0 vs cortical MAD | partial-volume | yes (TE, field, SWI phase-mask power, magnitude vs processed SWI, coil bias) | — | cortex mask from labels |
| Ridge percentile 85/80 within shell | — | — | — | — |
| `MIN_MM3`, `MIN_SLICES`, `MIN_ELONG`, `MIN_EXTENT_MM` | yes | blooming (TE/field) | `n_slices` = axis 2 | — |
| Merge gap = 2 × `vox[2]` | yes | — | assumes axis 2 = slice | — |
| `longs` = `n_sl × vox[2]` >40 mm "top-to-bottom" | yes | — | assumes axis 2 = S-I | — |
| Hessian σ 0.8/1.6 mm | anisotropic voxels | blooming | — | — |
| `pial_sd`, `cortex_frac`, `bank_frac`, `tram_frac`, E2, E3, E7 | yes (sub-voxel thresholds) | — | — | **yes, directly** |
| Mirror plane | — | — | head tilt | hemisphere labels |
| `edge`, rim | `vox[0]` only for the 2-D erosion | — | 2-D in array slices | — |
| ICH auto-finder | — | WM z | — | WM labels |
| FLAIR z | — | FLAIR protocol | — | FLAIR↔SWI rigid |

There is no machine-readable record of registration quality. `prep_anat.sh` prints the bbregister
`min cost` to the terminal only, and the mri_coreg route reports nothing.

### 1.5 Bugs and conceptual problems found

Severity: **H** = can change which candidates are shown or excluded on real data; **M** = biases a feature or
a validation number; **L** = fragility or cosmetic.

| ID | Sev | Where | Problem | Evidence |
|---|---|---|---|---|
| B1 | M | `detect_css.py:386, 480` | `extent_mm` uses `pad` slice widths minus 6. When `s.start−3` is clipped to 0, the width is underestimated by up to 3 voxels (7.5 mm on 2.5 mm slices). A real lesion on the first slices of a slab can be excluded as a "speck" (E4). Partial-brain slabs are exactly where this happens | code reading; certain |
| B2 | H | `detect_css.py:94, 113, 117, 323` | `CSF = [24]`. In recon-all `aparc+aseg`, sulcal CSF is label **0** inside the brain mask and label 24 is a small basal region (to be confirmed on P006: `np.unique(seg[brain & ~tissue])`). On the `--recon` route (a) `surface`/`surface_contact` miss sulcal CSF, (b) the FLAIR CSF reference is basal-cistern CSF or empty, and **FLAIR features silently switch off** when empty (`FL = None`). `tissue` itself is correct on both routes. Every test uses SynthSeg-style labels with 24 = extracerebral CSF, so the suite cannot catch this | code reading + FreeSurfer label convention; needs one-line check on P006 |
| B3 | M | `detect_css.py:220, 251` | 3-D Hessian is computed on `Ic` (0 outside the mask), not on `If` (mask filled with the cortical median, which Sato uses). The step at the mask edge gives large positive curvature that reaches 2σ ≈ 3.2 mm inside. Affects `tube_ratio`, `surface_alignment` and the vessel map near the convexity, especially with tight SynthStrip masks or thin extra-axial CSF | code reading |
| B4 | M | `detect_css.py:284-288` | `zmin3` is a 3×3×3 minimum of `z`, and `z` outside the brain is strongly negative (Ic=0). A mirror point that falls at or outside the mask edge (asymmetric head, tilt, convexity) reads as "dark". This could inflate `mirror_dark_frac` for the most superficial candidates. **Not reproduced on the phantom** (perfectly symmetric: 0 % of mirror points left the brain), so it is a hypothesis to check on P006/P011 | code reading; phantom negative |
| B5 | H | `detect_css.py:382-399` | Generation gates G1-G4 drop candidates without a record. `stress_test.py` reports these as "not detected", mixing "never dark enough" with "dropped by a gate". Violates the project's own rule that exclusions are listed with a reason | code reading |
| B6 | M | `detect_css.py:578-598` (E8) | Artifact-zone exclusion runs on the SWI-only route, unlike E2/E3 | code reading |
| B7 | M | `detect_css.py:179-203` | Merging is done before exclusion and features use the median over the merged object. A true cSS piece merged with an adjacent vein can be excluded as "vein" as a whole. `n_pieces` is computed but not written to the CSV | code reading (v4.9.1 saw the chain version of this on P006) |
| B8 | M | whole file | "Slices" means array axis 2. `mri_convert` preserves the stored orientation, so a sagittally or coronally stored SWI runs Sato, `n_slices`, merge gap and `longs` in the wrong plane. No check or warning | code reading |
| B9 | L | `detect_css.py:112, 485` | The global brain-edge mask `edge` is overwritten by the per-candidate scalar `edge` inside the loop. Harmless now (the mask is not used after line 146) but fragile. Also `shell |= edge` adds only voxels later removed by `~rim`; it does raise the ridge percentile | code reading |
| B10 | L | `detect_css.py:393` | `'n_pieces' not in dir()` at module scope as a flag. Works, fragile | code reading |
| B11 | M | `score_css.py:72-73`, detect FLAIR, `review_css.py` phase | `a2009s_swispace`, `flair_swispace` and phase files are not checked for staleness against the SWI (only `t1seg_swispace` is). Re-importing a SWI under the same ID keeps old Destrieux/FLAIR | code reading |
| B12 | L | `score_css.py:242` | Surface area = voxel count × V^(2/3). For anisotropic voxels (0.5×0.5×2 mm) the face area depends on orientation and the staircase overestimates area. Fine within one protocol, not across protocols (the very use case it was added for) | code reading |
| B13 | M | `feature_report.py:21-23` | With several experts, only the alphabetically first expert file is used. No consensus/adjudication | code reading |
| B14 | M | `feature_report.py`, `stress_test.py` | AUCs pool candidates across patients (no grouping), no confidence intervals | code reading |
| L1 | H (validation) | `rule_test.py` → v4.12, `CLAUDE.md` | **In-sample threshold selection.** E7 thresholds were chosen on P006 known foci vs P011 FPs and then reported as "P006 6/6 kept, P011 26→8" on the same data. That result is training performance, not validation | `CLAUDE.md` v4.12 |
| L2 | H (validation) | stress test on P006 | P006 is cSS-positive. Its real cSS count as "non-lesions" in the synthetic stress test, and its veins/cSS were used to design the synthetic generator | `CLAUDE.md` P006 stress |
| L3 | H (labels) | P006 "known" foci | The 6-7 known foci are v3 accepts by a non-expert reader who, blinded, later accepted 1/100 and called the old foci Normal/Vein. They are not a reference standard | `CLAUDE.md` |
| L4 | M | `match_calls` | Fingerprint tolerance 2 voxels / 15 % volume. A detector change that merges pieces drops the call (baseline: 6 of 8 calls dropped after a z −2.3 re-run). Safe direction, but it silently shrinks labelled data across versions | baseline log |

Checked and **not** a bug: anisotropic EDT (`sampling=vox` everywhere), normal-vector units in the
tram walk (`n / vox` converts mm to voxel steps), Hessian unit scaling (`/(vox[a]·vox[b]) · σ²`),
hemisphere assignment (majority of own voxels on the distance-based map), call re-attachment by
fingerprint.

### 1.6 Documentation vs code

| # | Document says | Code does |
|---|---|---|
| D1 | `PIPELINE.md` §7.6: score ×0.8 when mirror >0.6 | Not applied; mirror is recorded only |
| D2 | `PIPELINE.md` §5a table and `CLAUDE.md` v4.8: microbleed-like candidates are excluded | Flag only (`cmb_like`) since v4.9.2 |
| D3 | `PIPELINE.md` §5a and `CLAUDE.md` v4.8: artifact zone in the "lowest 30 %" | `rel_height < 0.4` (v4.9 text says 40 %) |
| D4 | `detect_css.py` docstring: "v4.11"; definition table lists `mirror_dark_frac` and `vein_tree_mm` as criteria | Package is v4.13; both are recorded only |
| D5 | `PIPELINE.md` §6 step 2 expects `vein_tree_mm` and `mirror_dark_frac` to be useful (<0.5) | Both demoted to recorded-only because they fail on real data |
| D6 | `PIPELINE.md` §7.5: FLAIR = median CSF z | Also cortical FLAIR z (`flair_ctx_*`) |
| D7 | `CLAUDE.md`: "on both phantoms every definition feature separates cSS from veins" | Baseline PH3: 2 of 7 mimics are still shown (central sulcal vein 1, surface vein 6). `check_v4.py` prints this but does not assert it |
| D8 | `CLAUDE.md` conventions: numpy/scipy/nibabel only | Also scikit-image, pandas, matplotlib |
| D9 | Brief for this task: "P006 ≈ 66 candidates after v4.12" | `CLAUDE.md`: 66 after v4.11, **59** after v4.12 |
| D10 | Brief: "SWI magnitude is the detection image" | `PIPELINE.md` §2 accepts "processed SWI image, or the SWI magnitude". These have different contrast (phase-mask weighting darkens veins); mixing them is a protocol confound that is not recorded per case |
| D11 | Brief cites Boston v2.0 as DOI `10.1016/S1474-4422(22)00220-9` | Search results give `10.1016/S1474-4422(22)00208-3`; check before citing |

### 1.7 Which tests mean something

Meaningful (deterministic, assert a specific behaviour): fingerprint re-attachment; ICH sulcal rule;
Boston v2.0 count; partial-slab coverage switch (PH5); cerebellar line excluded with cSS kept (PH4);
sulcal 3/4 score on PH3; agreement ICC on identical sessions.

Too easy or not asserting:
- **Contrast far above threshold.** Phantom cSS have `darkness_z` ≈ 7-10 (PH4 stress median 7.8); the detector threshold is 2.5. Nothing exercises faint lesions.
- **Geometry.** Ellipsoid brain, constant tissue intensities, Gaussian noise, no bias field, perfectly symmetric (mirror can never fail), straight isolated veins, uniform 2.5 mm cortex, open sulci.
- **v3 phantom.** Veins are flat planes in the banks; pass threshold is 4/8.
- **Stress tests print only.** Baseline PH4 with T1 anatomy: **5/8 found, AUC 0.40**, and the suite still passes.
- **`check_v4.py`** requires `found ≥ n−1` and AUC ≥0.8 on 5 cSS vs 2 mimics; surface veins being shown is not a failure.
- **Labels.** All phantoms use SynthSeg-style labels; recon-all `aparc+aseg` (B2) and Destrieux on the real grid are never tested.
- **No per-lesion reason assertion** except via `check_known` on PH4.

---

## 2. Literature evidence table

**Verification limit, stated plainly.** This environment's network policy blocks PubMed Central,
journal sites, DOI resolvers, arXiv and Europe PMC (HTTP 403 at the proxy). Only the web-search tool
works. Every row below was checked against search-result abstracts and snippets, not against full
text. Your brief asked for no reliance on snippets. I could not meet that, so the column "verified"
says what was actually seen, and rows marked "prior knowledge" are from my training and were not
re-checked this session. To lift the limit, allow those hosts under Network access in the cloud
environment settings (https://code.claude.com/docs/en/cloud-environments#network-access).

| Source | What it does for cSS | Sequences | Separates cSS from veins on real data? | Design / n | Reported metrics | Limits / implementable here? | Verified |
|---|---|---|---|---|---|---|---|
| Charidimou et al., *Neurology* 2017;89:2128. PubMed 29070669 | Visual multifocality score 0-4 | T2*-GRE | Human reader | Prospective, 313 CAA | κ 0.87 (score); annual ICH recurrence by score 0-4: 5, 6.5, 13.5, 16.2, 26.9 % | Defines the target; already implemented in `score_css.py` | abstract/snippet |
| Charidimou, *AJNR* 2016;37:E43. doi:10.3174/ajnr.A4748 | Minimum rating/reporting standards (letter) | T2*/SWI | Human reader; lists mimics | Letter | — | Criteria list not visible in snippets; `PIPELINE.md` quotes it | abstract/snippet (criteria not seen) |
| Charidimou et al., *Brain* 2015;138:2126 | Review: detection and significance | T2*/SWI | Review of mimics | Review | — | Supratentorial distribution; sparing of brainstem/cerebellum | search hit only |
| Charidimou et al., Boston criteria v2.0, *Lancet Neurol* 2022;21:714. doi:10.1016/S1474-4422(22)00208-3 (see D11) | cSS counted by gyri as haemorrhagic lesions | T2*/SWI | Human reader | Derivation + temporal + geographic validation vs histopathology | Probable CAA: sens 64.5 %, spec 95.0 % (autopsy subgroup); AUC 0.797 / 0.910 / 0.808. No cSS-alone accuracy found | Implemented as `boston2_css_lesions` | abstract/snippet |
| van Harten et al., *NeuroImage Clin* 2023;38:103447. doi:10.1016/j.nicl.2023.103447 | **Semi-automatic quantification** (not detection) | SWI | No: vesselness also catches veins; a human removes FPs | 20 CAA patients with cSS | Inter-observer Pearson 0.991; intra-observer ICC 0.995. 2-D Frangi, σ = 1 voxel | Seeds, region growing and Dice 0.75 are quoted in `PIPELINE.md` but were **not visible** in snippets this session | abstract/snippet |
| Cao et al., arXiv:2610.01542 (Oct 2026, MICCAI SASHIMI workshop) | **Automated cSS (and CMB) segmentation**, trained on synthetic data | MRI, synthesised from parcellations | Vessel mimics synthesised for CMB; cSS uses a hypointensity constraint | Test: 10 cSS, 13 CMB cases with manual delineations | cSS: AUPRC 0.284 vs Frangi 0.083; AUROC 0.907 vs 0.731 | Preprint, not peer reviewed, n=10. AUPRC 0.28 means most predicted voxels are false. A research comparator, not a validated detector. Note `make_synthetic.py` cites "Cao et al. 2026" for vein decoys; the snippets confirm vessel mimics only for CMB | abstract/snippet |
| Weidauer, Neuhaus, Hattingen, *Clin Neuroradiol* 2023;33:293. doi:10.1007/s00062-022-01231-5 | Review: aetiology, imaging | SWI/T2*, FLAIR | Review | Review | — | Context for mimics, infratentorial vs cortical SS, CAA-ri/ARIA-H | abstract/snippet |
| SWI vs GRE-T2* in advanced CAA, *Rev Neurol* 2023 (ScienceDirect S0035378723011384) | Visual rating, sequence comparison | SWI vs GRE | Human | 54 patients (38 vs 36 cSS-positive) | 70.4 % vs 66.7 %, p = 0.5; SWI showed more widespread cSS | Supports recording sequence and not pooling | abstract/snippet |
| 3D FLAIR / DIR for cSS, *Eur Radiol* 2021. doi:10.1007/s00330-021-07751-x | Sulcal hyperintensity related to cSS | FLAIR, DIR, SWI | — | Cognitive-clinic cohort | Less sensitive than SWI | Chronic cSS can show sulcal FLAIR signal; `flair_bright` must stay a hint, not an exclusion | abstract/snippet |
| Acute cSAH and cSS in CAA, *JNNP* 2018;89:397. doi:10.1136/jnnp-2017-316368 | Acute cSAH imaging and evolution | FLAIR, GRE | — | Probable CAA cohort | Acute cSAH = sulcal FLAIR hyperintensity + GRE hypointensity in all; evolved into cSS in all; TFNE 76 % vs 34 % | Supports FLAIR as acute-cSAH flag | abstract/snippet (authors not checked) |
| UCL ISMRM 2019 abstract #2933 (BOCAA cohort) | SWI non-local phase effects | SWI, QSM | — | Conference abstract | Qualitative: SWI broadens/duplicates CMB and deforms superficial siderosis; QSM delineates better | No quantitative validation; QSM not used (site decision) | abstract/snippet |
| Al-Masni et al., *NeuroImage Clin* 2020 (CMB, not cSS) | Two-stage CMB detection | SWI + phase | Veins, calcium as mimics | Multi-subject | 93.6 % stage-1 sensitivity; 1.42 FP/subject after stage 2 | Shows the candidate→FP-reduction design and FP/subject reporting; CMB, not cSS | abstract/snippet |
| QSM-based CMB detection, *JMRI* doi:10.1002/jmri.29198 | Two-stage CMB vs mimics | QSM | Vessel mimics | — | Sens 88.9 %, 2.87 FP/subject | CMB, not cSS | abstract/snippet |
| Riley et al., *BMJ* 2020;368:m441 | Sample size for prediction models | — | — | Methods | Replaces the 10-EPV rule with shrinkage/optimism criteria (`pmsampsize`) | Use to size the logistic model before fitting | abstract/snippet |
| Subject-level leakage: arXiv 2309.00350 (brain MRI) and others | Leakage from slice/image-level splits | — | — | Mostly preprints | Inflation varies widely | Supports patient-grouped CV | abstract/snippet |
| Frangi et al., MICCAI 1998; Sato et al., *Med Image Anal* 1998;2:143 | Hessian vesselness | — | Known to respond to any dark ridge, incl. veins | Methods | — | Already used | prior knowledge |
| Haacke et al., *AJNR* 2009;30:19 (SWI technical review), doi:10.3174/ajnr.A1400 | Phase handedness, calcium vs blood | SWI phase | Phase separates dia- from paramagnetic, not haemosiderin from deoxy-Hb | Review | — | Supports the existing reading-aid use of phase | prior knowledge |
| Greve & Fischl, *NeuroImage* 2009;48:63 (bbregister) | Boundary-based registration cost | T1 + EPI/T2* | — | Methods | — | Record `mincost` per case as QC | prior knowledge |
| FreeSurfer in AD/iNPH with manual edits, *Front Neurosci* 2024. doi:10.3389/fnins.2024.1366029 | Segmentation errors in atrophy | T1 | — | 19 iNPH, 28 AD, 30 controls | Edits changed volumes more in iNPH | No pial-surface error rates for elderly found | abstract/snippet |

**Searched for and not found:** a peer-reviewed, real-world-validated, automated cSS detector; any
cSS method in cortical-surface or geodesic coordinates; any study measuring vesselness false positives
specifically for cortical veins along the pial surface in cSS; a multi-vendor or 1.5 T vs 3 T cSS
detection study; pial-surface accuracy near the tentorium in elderly brains. The closest automated
work is the Cao 2026 preprint (n=10 cSS, AUPRC 0.28). On current evidence **no robust real-world
automated cSS detector exists**; human review remains the published standard.

---

## 3. Ranked false-positive reduction plan

Ranking = expected value ÷ risk to sensitivity. "Hard exclusion" is reserved for mimics that are
specific and validated; everything new starts as a recorded feature, reviewer flag or tier.

### 3.1 Safe immediate changes

**S1. Per-case QC JSON and a case confidence state**
```
Recommendation: write review/ID_qc.json and a case_status of ok / low_confidence_review /
  insufficient_quality; append qc_status to every candidate row.
Why it should reduce false positives: it does not remove any candidate. It stops anatomy-based
  rules (E2, E3, E7, E8) from making confident decisions on bad anatomy, and tells the reader
  which cases' candidate lists are untrustworthy.
Expected sensitivity risk: none (no candidate is removed; a low-confidence case gets fewer
  anatomy-based exclusions, i.e. more candidates).
Evidence: CLAUDE.md P010 (SWI-only rules backfired); van Harten 2023 and the SWI vs GRE study
  (protocol dependence); Greve & Fischl 2009 (bbregister cost).
Data needed: none to build. Thresholds for "low confidence" need ~10 labelled cases.
Implementation location: new scripts/css_qc.py (called from run_css.sh after align_seg.py, and
  importable by detect_css.py); prep_anat.sh writes the bbregister/mri_coreg cost to work/ID_reg.json.
Validation test: phantom with label map shifted 0/1/2/3 mm → registration-edge score degrades
  monotonically and status flips at the configured shift; coronal/sagittal-stored copy → orientation
  warning.
Decision rule: QC failure (case-level) / reviewer flag.
```

**S2. Bug fixes B1, B3, B5, B6, B9, B10, B11 (+ B2 after a one-line check on P006)**
```
Recommendation: fix extent clipping (B1), Hessian input (B3: compute on If), record generation
  drops with reasons (B5: work/ID_dropped.csv, not shown to the reader), gate E8 on T1 like E2/E3
  (B6), rename the scalar edge (B9), carry n_pieces explicitly (B10), mtime-check a2009s/FLAIR/phase
  (B11). B2: define CSF as brain ∧ ¬tissue (covers label 0 inside the mask and label 24) after
  confirming the label distribution on P006.
Why it should reduce false positives: B3 removes mask-edge Hessian artefacts that make convexity
  noise look like sheets/tubes; B2 restores FLAIR features on the recon-all route. B1, B5, B6 mainly
  protect sensitivity and make misses explainable.
Expected sensitivity risk: low; B3 and B2 change feature values, so known foci must be re-checked.
Evidence: code audit, baseline tests.
Data needed: P006 and P011 re-runs with check_known.py.
Implementation location: detect_css.py lines 94-117, 220, 251, 386-399, 480, 597; score_css.py 72.
Validation test: unit tests per bug (section 7); run_tests.sh; check_known P006 6/6; P011 count
  reported before/after.
Decision rule: bug fix (no new rule).
```

**S3. Two-tier review list instead of more deletions**
```
Recommendation: add review_tier (1 = likely cSS, 2 = low-confidence, not safe to exclude). Tier 2 =
  any candidate with a mimic flag (artifact_zone, vein_like, near sinus, cmb_like, low qc) or
  css_evidence below a reader-chosen cut. Report tier-1 and tier-2 counts per scan. Keep the
  excluded file unchanged.
Why it should reduce false positives: it lowers the effective review burden (readers can read
  tier 1 first and see how many tier-2 remain) without hiding anything.
Expected sensitivity risk: none at the candidate level; there is a reader-attention risk if tier 2
  is skipped, which is why it is reported, never dropped.
Evidence: two-stage CMB literature (candidate stage + FP-reduction stage, FP/subject reporting).
Data needed: tier rules can be set now as transparent rules; the cut is tuned only after labels.
Implementation location: detect_css.py after exclusion(); review_css.py sorts tier 1 first.
Validation test: phantom → every truth lesion in tier 1 or 2; P006 known foci tier recorded.
Decision rule: soft score / reviewer flag.
```

**S4. `rule_version` and `exclusion_confidence` on excluded rows**
```
Recommendation: append rule_version (e.g. "E7@v4.12") and exclusion_confidence (high for E1/E4 on
  T1, medium for E5/E6/E7, low for any rule fired on a low-QC case) to ID_excluded.csv.
Why it should reduce false positives: none directly; it makes low-confidence exclusions easy to
  audit and re-include.
Expected sensitivity risk: none.
Decision rule: reviewer flag.
```

### 3.2 Changes that need expert-labelled real data before any threshold

| ID | Change | Why not now |
|---|---|---|
| D-1 | Orbitofrontal / temporal-pole artifact rule above `rel_height` 0.4 (5 of 8 remaining P011 FPs) | One negative patient; artifact extent depends on sinus anatomy, TE and field. Needs ≥10 negatives with expert "artifact" labels and positives in those regions |
| D-2 | Turn E3/E8 (skull base) from hard exclusion into tier 2 | Increases candidates; needs labelled inferior temporal/occipital cSS to judge the trade-off. This changes the phenotype the detector can find, so it is a scientific decision for the study team |
| D-3 | Any threshold on `mirror_dark_frac`, `vessel_run_mm`, `vessel_ext_mm`, new vein-tracking features | Real-data AUCs so far are ~0.3-0.5 or synthetic-artifact-driven |
| D-4 | Regularised logistic ranking (replaces `score_v4` ordering only) | Needs ~20-30 positive and 20-30 negative patients (Riley 2020 criteria via `pmsampsize`, events counted per **patient**, not per candidate) |
| D-5 | Re-deriving E7 with patient-grouped CV | Current thresholds are in-sample (L1) |

### 3.3 Research prototypes (run beside the detector, write extra columns, never remove)

| ID | Idea (brief section) | Assessment |
|---|---|---|
| R-A | **Surface/geodesic representation** (A). Map candidate voxels to `lh/rh.pial` and `white` vertices (nibabel `read_geometry`, SWI→T1 via the LTA); graph on mesh edges; geodesic length (Dijkstra, `scipy.sparse.csgraph`), geodesic curvature of the vertex-path centreline, covered area, normal-distance distribution (signed distance along vertex normals), fraction explained by a surface-parallel band, banks/crown/fundus coverage from `sulc`/`curv`, abrupt end vs continuation | Highest-value prototype: models the definition directly and avoids Euclidean distances across folds. Needs recon-all (Route B only) and good registration (S1). No published cSS precedent; must be shown to add over voxel features |
| R-B | **Local bounded vein tracking** (B). Seed at candidate endpoints; minimal-path / tracking on a cost = f(vesselness, z, direction change), radius 15-25 mm; record continuation per end, branch count, diameter CV along the path, orthogonal-plane "dot" test (2-D cross-section roundness), tube-vs-sheet along the path, distance to a sinus mask (SynthSeg has no sinus label: needs a sinus atlas or a dark-tube-at-midline/transverse-region heuristic) | Replaces the meaning of `vein_tree_mm`. Main risk: cSS next to a vein inherits the vein's continuation; keep per-end values and the gap between candidate and vessel |
| R-C | **Geodesic tram-track** (C). Opposite bank found on the pial mesh across the sulcal fundus (vertex pairs across a sulcus from `sulc` and normals), not by one voxel ray | Must report `single_bank` as a valid outcome, never as evidence against cSS |
| R-D | **Multi-sequence mimic flags** (D). Phase: local filtered-phase sign vs a reference vein on the same image (vendor handedness from JSON `Manufacturer`); FLAIR: acute cSAH; DWI/ADC: cortical restriction; T1: cortical hyperintensity (laminar necrosis) | Phase separates diamagnetic calcium from paramagnetic material only; it does **not** separate haemosiderin from deoxy-Hb. Each flag is reviewer evidence |
| R-E | **Synthetic-trained segmentation model** (Cao 2026 style) as a comparator | Preprint; interesting for candidate generation, not a replacement |

---

## 4. Data schema for expert labels and hard negatives

All files under `$CSS_BASE`, study IDs only, never committed.

`labels/candidate_labels.csv` (one row per reader × candidate × session)

| column | type | notes |
|---|---|---|
| study_id | str | P006 |
| reader_id | str | expertA |
| session | int | 1, 2 (intra-rater) |
| detector_version | str | v4.13 + git hash |
| cand_uid | str | `study_id:vox_i:vox_j:vox_k:volume_mm3` fingerprint (stable across re-runs within tolerance) |
| list | enum | candidate / excluded / dropped (so excluded and dropped items can be labelled too) |
| label | enum | cSS, vein_surface, vein_sulcal, vein_transcortical, thrombosed_vein, acute_cSAH, microbleed, ICH_related, calcification, artifact_airbone, artifact_motion, normal_iron_cortex, laminar_necrosis, haemorrhagic_infarct, infratentorial_SS, uncertain |
| confidence | int 1-5 | |
| sulcus / gyrus | str | Destrieux name if known |
| blinded | bool | reader saw no scores/features |
| timestamp | ISO 8601 | |

`labels/missed_lesions.nii.gz` + `labels/missed_lesions.csv`: reader-drawn cSS that **no candidate
covers**. Without this, candidate-level labels can only measure precision, never sensitivity.

`labels/consensus.csv`: adjudicated label per `cand_uid` and per missed lesion, with
`adjudication_method` (agreement / discussion / third reader).

`labels/study_meta.csv` (one row per study): cSS_status (pos/neg/uncertain, patient-level reference),
reference multifocality 0-4 per hemisphere, vendor, field_T, TE_ms, TR_ms, sequence (SWI-magnitude /
SWI-processed / GRE), acq_dim (2D/3D), voxel_mm, slice_mm, coverage (full / slab / no skull base),
motion_grade 0-3, T1_available, recon_all, FLAIR, phase, DWI, qc_status (from S1), split (dev/test,
assigned once, see 5).

Hard-negative set = all consensus labels other than cSS and uncertain, with `vein_surface`,
`vein_sulcal`, sinus-adjacent veins, `artifact_airbone` and `normal_iron_cortex` as named strata so
they can be over-sampled in ranking evaluation (never in the test-set metrics).

---

## 5. Leakage-safe validation protocol

1. **Split by patient, once, before looking.** Assign each new study to `dev` or `test` at enrolment
   (stratified by cSS status and protocol). P006, P010, P011 and every study used to design rules so
   far are `dev` forever.
2. **Candidate generation is label-blind.** Detector parameters are frozen per version; no label,
   reader call or post-review quantity (accepted mask, grown volume, sulcal assignment) is an input.
3. **Model selection inside dev only**: leave-one-patient-out or grouped K-fold (K=5) with nested CV
   for regularisation strength. All candidates of a patient are in the same fold. Feature
   standardisation fitted inside each training fold.
4. **Calibration** inside the nested loop (Platt or isotonic with care at small n); report calibration
   slope/intercept and Brier score, with bootstrap optimism correction.
5. **Test once.** The frozen pipeline (detector + ranking + tiers) runs on `test` once per release.
6. **Metrics**, all with patient-bootstrap 95 % CIs:
   - candidate sensitivity (against consensus cSS plus missed-lesion annotations) and FP per scan;
   - FP per **negative** scan;
   - sensitivity at review burden top-5/10/20 and tier 1 only;
   - FROC (sensitivity vs mean FP per scan) and PR-AUC;
   - patient-level sensitivity/specificity for "any cSS";
   - multifocality score agreement (weighted κ) and surface-area ICC vs expert;
   - stratified by vendor, field, sequence, slice thickness, T1 route, coverage, qc_status.
7. **Promotion rule for a new hard exclusion**: on `test`, it removes **no** adjudicated cSS and
   reduces FP per negative scan with a CI that excludes zero. Otherwise it stays a flag/tier.
8. **Report AUC only alongside FP per scan.** A ranking AUC of 0.9 with 30 candidates per negative
   scan is not usable.

---

## 6. Staged implementation plan

Each stage is one commit, runs `bash tests/run_tests.sh` + the new tests, and on real data
`check_known.py P006` and a P011 count. Rollback for every stage = `git revert <commit>`; no stage
deletes or renames an existing column, file or rule.

### Stage 1 (smallest safe change, proposed for approval first): QC JSON + bookkeeping, no candidate changes

- **Files**: new `scripts/css_qc.py`; `scripts/detect_css.py` (end of file: output section);
  `scripts/run_css.sh` (call after `align_seg.py`); `scripts/prep_anat.sh` (save registration cost).
- **Pseudocode**
  ```
  qc = {}
  qc.protocol  = read raw/ID/SWI/*.json (vendor, field, TE, TR, slice, ImageType -> magnitude/processed/GRE)
  qc.voxel_mm, qc.orientation = aff2axcodes(swi.affine); qc.slice_axis = argmax |affine col| on S-I
  qc.orientation_ok = (slice_axis == 2)
  qc.coverage = FULL_BOTTOM logic + top-of-head check (brain area at top slice)
  qc.label_coverage = labelled ∩ brain / brain
  qc.reg_cost = work/ID_reg.json (bbregister mincost or mri_coreg final cost) or null
  qc.edge_agreement = median |pial_sd| at voxels of strong SWI gradient magnitude near cortex
  qc.edge_agreement_by_region = same per lobe x hemisphere
  qc.motion_proxy = ratio of background (outside-head) signal in the raw SWI vs brain, if raw kept
  qc.anatomy = T1 recon | T1 synthseg | SWI-only
  status = insufficient_quality if orientation_ok is false or label_coverage < x
           low_confidence_review if anatomy is SWI-only or edge_agreement > y or reg_cost > z
           else ok
  ```
  x, y, z start as placeholders reported but **not used to change exclusions** in Stage 1.
- **New outputs**: `review/ID_qc.json`; appended candidate columns `qc_status`, `rule_version`,
  `n_pieces`; appended excluded columns `rule_version`, `exclusion_confidence`.
- **Also in Stage 1** (output-neutral or sensitivity-protecting fixes): B1 extent clip, B5 dropped-candidate log `work/ID_dropped.csv`, B9/B10 renames.
- **Tests**: phantom label shift 0/1/2/3 mm → `edge_agreement` increases monotonically; transposed
  (sagittal-stored) phantom → `orientation_ok=false`; candidate on slice 0 of a slab → `extent_mm`
  equals the in-volume value; generation drops appear in `ID_dropped.csv` with reason.
- **Acceptance**: all existing tests pass; existing CSV columns identical in name, order and value
  except `extent_mm` for edge-clipped candidates; P006 6/6 kept; P011 candidate count unchanged
  unless an edge-clipped speck is restored (reported).
- **Decision type**: QC output + reviewer flags. No exclusion behaviour changes.

### Stage 2: feature-changing bug fixes (B2, B3, B6, B11)

- **Files/regions**: `detect_css.py` 94-117 (CSF definition), 215-260 (Hessian input `If`), 597 (E8
  gated on T1); `score_css.py` 72 and `detect_css.py` 319 (staleness checks).
- **Pre-check**: on P006 recon-all labels, print label counts in `brain ∧ ¬tissue`.
- **Tests**: phantom variant with aparc+aseg convention (sulcal CSF = 0, label 24 only basal) →
  `surface_contact` and FLAIR features equal the SynthSeg-convention values within tolerance;
  Hessian at mask edge: synthetic flat image with a tight mask → no tube/sheet response within 3 mm.
- **Acceptance**: P006 6/6 known kept; changes in P006/P011 candidate counts reported per rule; no
  phantom lesion lost.
- **Rollback**: revert; outputs are regenerated by re-running detect.
- **Decision type**: bug fix.

### Stage 3: two-tier review list (S3)

- **Files**: `detect_css.py` after `exclusion()`; `review_css.py` ordering; `export_review.py`
  unchanged (columns appended after the stable block).
- **Pseudocode**: `tier = 2 if (qc_status != ok) or artifact_zone or vein_like or cmb_like or
  near_sinus or css_evidence < CSS_TIER_CUT else 1` with `CSS_TIER_CUT` an env var default that puts
  every phantom lesion in tier 1.
- **Columns**: `review_tier`. Summary line: `tier1 N, tier2 M`.
- **Acceptance**: no phantom/synthetic lesion in tier 2 at default; tier counts reported for P006/P011.
- **Decision type**: soft score / reviewer flag.

### Stage 4: recorded-only vein-tracking features (R-B) and phantom mimics (section 7)

- **Files**: new `scripts/css_vessel.py` (functions), called from `detect_css.py` feature loop.
- **Columns**: `local_vein_cont_end1_mm`, `local_vein_cont_end2_mm`, `local_vein_branch_count`,
  `local_vein_diam_cv`, `orthogonal_dot_frac`, `sinus_dist_mm`, `sinus_direction_score`.
- **Thresholds**: none. Recorded only.
- **Acceptance**: on new phantom surface-following veins, AUC vs cSS reported; no change to candidates.

### Stage 5: surface/geodesic prototype (R-A, R-C), recon-all route only

- **Files**: new `scripts/css_surface.py`; reads `subjects/ID/surf/{lh,rh}.{pial,white,sulc,curv}`
  and `work/ID_swi2t1.lta`.
- **Columns**: `geo_length_mm`, `geo_area_mm2`, `geo_curvature`, `normal_dist_med_mm`,
  `normal_dist_iqr_mm`, `band_explained_frac`, `banks_covered` (0/1/2), `crosses_fundus`,
  `geo_tram`, `single_bank`.
- **Acceptance**: recorded only until stage 6.

### Stage 6: patient-grouped feature analysis and logistic benchmark (needs ≥20+20 labelled patients)

- **Files**: new `scripts/fit_rank.py` (sklearn would be a new dependency: ask first) or a small
  numpy L2-logistic implementation.
- **Output**: `results/rank_cv.json` with grouped-CV FROC, PR-AUC, calibration, per-feature
  coefficients with bootstrap CIs; comparison: `score_v4`, logistic, tier-only.
- **Promotion**: only by the rule in section 5 #7.

---

## 7. Test-suite gaps and additions

Each new synthetic object gets an id and the test asserts its **fate** (shown / excluded with reason
X / dropped with reason Y / not detected), not just aggregate sensitivity.

| New phantom content | Asserts |
|---|---|
| Cortical vein running **along** the pial surface for 30-50 mm, 0.5-1 mm outside, curved | Reported fate and feature values; becomes a regression target for stages 4-5 |
| Vein that branches and joins a midline "sinus" tube | `local_vein_*`, `sinus_*` features |
| Asymmetric veins; bilateral symmetric cSS | Mirror cannot be used as evidence either way |
| Thin single-bank cSS (1 voxel) | Shown; `single_bank` |
| Curved cSS from crown into fundus | Not excluded by E6 (`axis_normal`) |
| cSS at pial depth −1.5, −0.5, +0.5 mm | Shown at all three |
| Same lesions at 0.5×0.5×1, 0.9 isotropic, 0.5×0.5×2.5 mm, with blooming 0/1/2 voxels | Per-resolution fate table |
| Partial slab starting mid-brain, lesion on the first 2 slices | Not excluded as speck (B1) |
| Orbitofrontal/temporal-pole signal dropout field | Fate + `artifact_zone` |
| Phase phantom with calcium-sign and blood-sign blobs, both vendor handedness conventions | Phase flag uses a same-image vein reference |
| Acute cSAH: sulcal FLAIR-bright + SWI dark | `flair_bright = 1` |
| Thrombosed vein with adjacent cortical FLAIR/DWI hyperintensity | `flair_ctx_bright = 1` |
| Laminar necrosis (T1-bright cortex), haemorrhagic infarct (blob in cortex + DWI) | Flags only |
| Label map shifted 1, 2, 3 mm | QC status changes; anatomy rules downgraded |
| Atrophy: widened sulci (4-6 mm CSF) and thinner cortex | Tram-track still found; E7 does not fire on cSS |
| Motion ghosting (shifted copy at 10 %) and low SNR (noise ×2, ×3) | QC motion proxy; FP count reported |
| aparc+aseg label convention (sulcal CSF = 0) | B2 regression |
| Faint lesions: darkness near threshold (z ≈ −2.5 to −3.5), not 7-10 | Fate table at threshold |
| Make stress-test outcomes assertive | e.g. PH4 T1 host sensitivity ≥ 7/8 at depth 0.6 (currently 5/8 passes silently) |

The phantoms remain regression tests. None of them stands in for expert-labelled real cases.

---

## 8. What remains unvalidated

- The v4 ranking, `css_evidence`, `vein_evidence` and every hard-exclusion threshold were set by hand
  on phantoms, one cSS-positive patient with non-expert labels (P006), and two negatives (P010
  SWI-only, P011). None has been tested on held-out patients.
- E7 (off-cortex) thresholds are in-sample (L1).
- Sensitivity to faint real cSS, to cSS over the tentorium/skull base (excluded by design), and to
  protocols other than the few seen (vendor, 1.5 T, GRE, thick slices) is unknown.
- FLAIR features may be inactive on the recon-all route (B2), pending a check.
- No inter-rater study of the pipeline output exists yet.
- The tool is research software for candidate generation with mandatory human review. It is not
  clinically validated and must not be used for clinical decisions.
