# cSS grading pipeline v5.1.1

Three-layer detection: definition-based candidates, a mimic layer (phase, FLAIR, microbleed clusters) and a ranking
model trained on expert labels. **v5.0.0 (Stage A):** both new layers are in place but OFF by default
(`CSS_USE_MIMIC=1`; the ranker needs >= 3 expert-labelled patients and is applied to dev-split studies only);
with the defaults the output equals v4.19. v5.1.0 adds per-object fate accounting, `audit_case.py`, stale-input
checks and a fix to the opt-in microbleed rule (docs/CHANGELOG_v5.1.md). v5.1.1: per-piece evidence of merged objects (`CSS_MIXED_REVIEW=1`), valid-sample fractions, missing-input indicators, label schema. The pial-surface sampler (layer 2) is still to come.

Research software, not clinically validated. **How to run it, step by step: [RUNBOOK.md](RUNBOOK.md).**
Methods and rationale: [PIPELINE.md](PIPELINE.md). Audit and plan: docs/AUDIT_v4.13_and_FP_plan.md.
Technical context for Claude Code: CLAUDE.md.

## Install on the Mac
    cd ~/Downloads/cSS_pipeline && zsh install.sh      # backs up old scripts first
Open a new Terminal. If needed once: `conda install -n css -c conda-forge matplotlib scikit-image -y`
Check: `bash tests/run_tests.sh` and `bash tests/test_shell.sh` must both pass.
v5.0 switches (all off by default): `CSS_USE_MIMIC=1`, `CSS_QC_STRICT=1` - see RUNBOOK section 13.

## Per patient
    dcm2niix -z y -f "%s_%d" -o ~/css_project/raw/P007/SWI   "<SWI DICOM folder>"     # NOT the minIP
    dcm2niix -z y -f "%s_%d" -o ~/css_project/raw/P007/T1    "<3-D T1 DICOM folder>"
    dcm2niix -z y -f "%s_%d" -o ~/css_project/raw/P007/FLAIR "<FLAIR DICOM folder>"
    run_css.sh P007 <swi.nii.gz> --t1 <t1.nii.gz> --flair <flair.nii.gz> [--recon]   # strip, anatomy, detect
    python ~/css_project/scripts/review_css.py P007                  # review ALL candidates + score
    python ~/css_project/scripts/feature_report.py                   # features vs your calls (all cases)
`--recon` runs FreeSurfer recon-all (1-3 h) and enables scoring by sulci (Destrieux).

## Contents
scripts/   pipeline scripts (installed to ~/css_project/scripts)
tests/     phantoms (v3 + v4 definition phantom), run_tests.sh, test_shell.sh (no patient data needed)
tools/     build_workbook.py (Excel validation workbook)
RUNBOOK.md   the full procedure, DICOM to score, with every command
PIPELINE.md  how to read it, method details, validation plan
CLAUDE.md  project context for Claude Code
ALL_CODE.md every script in one file for reading
