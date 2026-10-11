# cSS grading pipeline v4

Research software. **Step-by-step guide: [PIPELINE.md](PIPELINE.md).** Technical context for
Claude Code: CLAUDE.md.

## Install on the Mac
    cd ~/Downloads/cSS_pipeline && zsh install.sh      # backs up old scripts first
Open a new Terminal. If needed once: `conda install -n css -c conda-forge matplotlib scikit-image -y`
Check: `bash tests/run_tests.sh` and `bash tests/test_shell.sh` must both pass.

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
PIPELINE.md  how to run it, how to read it, method details, validation plan
CLAUDE.md  project context for Claude Code
ALL_CODE.md every script in one file for reading
