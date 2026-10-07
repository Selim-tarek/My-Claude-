# cSS grading pipeline v3.2

Research software. See CLAUDE.md for the full technical context.

## Install on the Mac
    cd ~/Downloads/cSS_pipeline_v3.2 && zsh install.sh      # backs up old scripts first
Open a new Terminal. If needed once: `conda install -n css -c conda-forge matplotlib -y`

## Per patient
    dcm2niix -z y -f "%s_%d" -o ~/css_project/raw/P007/SWI "<SWI DICOM folder>"
    run_css.sh P007 ~/css_project/raw/P007/SWI/<file>.nii.gz       # strip, segment, detect
    python ~/css_project/scripts/review_css.py P007 --top 20         # interactive review + score
    python ~/css_project/scripts/feature_report.py                   # features vs your calls

## Contents
scripts/   all pipeline scripts (installed to ~/css_project/scripts)
tests/     phantom brain + run_tests.sh (no patient data needed)
tools/     build_workbook.py (Excel validation workbook)
CLAUDE.md  project context for Claude Code
ALL_CODE.md every script in one file for reading

## Review in Claude Code
    cd ~/Downloads/cSS_pipeline_v3.2
    claude
Claude Code reads CLAUDE.md automatically. Suggested first prompt:
"Read CLAUDE.md, run tests/run_tests.sh, then review scripts/detect_css.py and scripts/score_css.py
for bugs and for the v4 plan."
