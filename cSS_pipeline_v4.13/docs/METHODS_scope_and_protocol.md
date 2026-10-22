# Methods text (draft): phenotype scope and imaging protocol

Draft for the PI / study protocol. The choices marked **[decide]** are study decisions, not software
defaults; the paragraph must state whichever is chosen. Pipeline version: v4.18.

## Phenotype scope

The detector targets **supratentorial convexity cSS**: curvilinear SWI hypointensity lining the cortical
surface of the cerebral convexities and sulci, as in the Boston criteria v2.0 and the Charidimou 2017
multifocality scale. By design, it does not report dark structures in the regions below. They are
removed by hard exclusion rules; every removed object is kept in `review/ID_excluded.csv` with the rule id:

| Rule | Region removed | Why | Cost (cSS that cannot be detected) |
|---|---|---|---|
| E1 | ≥ 30 % within 5 mm of cerebellum / brainstem (tentorial interface) | cSS is supratentorial by definition; tentorial veins and sinuses | cSS on the inferior temporal / occipital surface lying on the tentorium |
| E2 | ≥ 30 % within 5 mm of brainstem, ventral DC, hippocampus, amygdala (T1 anatomy only) | circle of Willis, basal veins | mesial temporal cSS |
| E3 | lowest 20 % of the cerebrum height (T1 anatomy and full coverage only) | skull-base susceptibility artifact | basal frontal / temporal cSS |
| E8 | orbitofrontal / temporal-pole / inferior-temporal labels in the lowest 40 % (T1 only) | frontal-sinus / mastoid air artifact | orbitofrontal and anterior temporal cSS |

**[decide]** Either:

- **(a) Keep the scope** and write: "Automated detection was restricted to supratentorial convexity cSS
  above the skull base; cSS on the basal temporal, basal occipital and orbitofrontal surfaces was not
  assessed by the software and was rated by the expert readers only." In this case the expert must still
  read the basal surfaces, and those foci are reported separately, not as software misses. Or:
- **(b) Include basal cSS**: E1/E2/E3/E8 become review flags shown in a low-priority tier instead of
  exclusions. This raises reader burden (on P011 (v4.17.1, before re-registration), 18 of 44 excluded objects came from these rules) and
  needs the evaluation to report sensitivity for basal and convexity cSS separately.

Infratentorial (cerebellar / brainstem) superficial siderosis is a different disease and is excluded under
both options.

## Imaging protocol

- **Required:** SWI (or T2*-GRE, recorded separately) magnitude, skull-stripped by the pipeline;
  **3-D T1** in the same session. The T1 gives the cortical labels that the position-based rules and the
  sulcal scoring rely on. Without it (SWI-only route), the anatomy comes from SynthSeg on the SWI. That is
  too coarse for these features (P010: feature AUCs ~0.5), so position rules E2/E3/E8 are switched off and
  the case is `low_confidence_review`.
- **Recommended:** FreeSurfer recon-all on the T1 (Destrieux sulcal labels for the multifocality score;
  boundary-based registration, `bbregister`). The intensity-based fallback (`mri_coreg`) produced a 3.9 mm
  up-down label slip on P011.
- **Registration QC:** every case is checked in the axial, coronal and sagittal planes, plus the
  automated rigid-shift estimate (`label_shift_mm`); cases > 1.5 mm are re-registered before reading.
- **Orientation:** the SWI must be stored with axial slices (array axis 3 = superior-inferior); the
  detector stops otherwise.
- **Voxel size:** recorded per case (`review/ID_qc.json`). Slices > 2 mm are `low_confidence_review`.
  Ridge-filter scales are fixed in mm (0.6 / 1.2 / 1.8 mm, v4.18), but distance thresholds near 0.3–0.5 mm
  (rules E5, E7) are below the voxel size of thick-slice data. Therefore **all performance figures are
  reported stratified by in-plane voxel size and slice thickness**.

## Reference standard and reporting

- Expert reading with fine label codes (cSS, surface / sulcal / transcortical / sinus veins, artifacts,
  mimics). Missed cSS are outlined, and two readers reach consensus with adjudication (`consensus.py`).
- The development cases (P006, P010, P011, on which thresholds were set) are excluded from all reported
  performance. A test split is assigned at enrolment and is never used for tuning.
- Primary metrics are patient-grouped: sensitivity per lesion and per patient, false positives per
  cSS-negative scan, and the agreement of the multifocality score (weighted κ), all with patient-bootstrap
  confidence intervals. Pooled candidate-level AUC is not a primary metric.
- Research software; not validated for clinical use.
