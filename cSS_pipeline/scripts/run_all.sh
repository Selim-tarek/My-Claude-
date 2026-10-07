#!/bin/zsh
# Batch: run several patients without opening the viewer.
# usage: run_all.sh P006 P007 P008        (expects data/<ID>_swi.nii already imported)
B=${CSS_BASE:-$HOME/css_project}
for s in "$@"; do
  echo "======== $s"
  NOVIEW=1 $B/scripts/run_css.sh $s | head -1
done
