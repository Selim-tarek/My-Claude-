#!/bin/zsh
# One patient, start to finish (up to review).
# usage: run_css.sh SUBJECT_ID [path/to/swi.nii(.gz)] [--nostrip]
#   - with a path: imports the scan; skull-strips it with SynthStrip unless --nostrip
#     (use --nostrip only for scans that are already skull-stripped)
#   - set NOVIEW=1 to skip opening freeview (batch use)
set -e
S=$1; SRC=$2; STRIP=1
[[ "$3" == "--nostrip" ]] && STRIP=0
B=${CSS_BASE:-$HOME/css_project}; SC=$B/scripts
if [ -z "$S" ]; then
  echo "usage: run_css.sh SUBJECT_ID [path/to/swi.nii or .nii.gz] [--nostrip]"; exit 1
fi
mkdir -p $B/data $B/synthseg $B/work $B/review
if [ -n "$SRC" ]; then
  if [ $STRIP -eq 1 ]; then
    echo "[0/4] importing + skull-stripping (SynthStrip, CSF kept) $SRC"
    mri_convert "$SRC" $B/work/${S}_raw.nii.gz > /dev/null
    mri_synthstrip -i $B/work/${S}_raw.nii.gz -o $B/data/${S}_swi.nii \
      -m $B/work/${S}_brainmask.nii.gz > /dev/null
  else
    echo "[0/4] importing (no skull-strip) $SRC"
    mri_convert "$SRC" $B/data/${S}_swi.nii > /dev/null
  fi
fi
if [ ! -f $B/data/${S}_swi.nii ]; then
  echo "ERROR: $B/data/${S}_swi.nii not found - give the scan path as 2nd argument"; exit 1
fi
echo "[1/4] SynthSeg segmentation"
# re-run when missing OR when the SWI was re-imported after the last segmentation (v4 fix:
# v3 silently reused a stale segmentation for a different scan with the same ID)
if [ ! -f $B/synthseg/${S}_seg.nii.gz ] || [ $B/data/${S}_swi.nii -nt $B/synthseg/${S}_seg.nii.gz ]; then
  mri_synthseg --i $B/data/${S}_swi.nii --o $B/synthseg/${S}_seg.nii.gz \
    --parc --robust --threads 8 > /dev/null
else
  echo "      (already done, skipping)"
fi
echo "[2/4] aligning segmentation to SWI"
python $SC/align_seg.py $S > /dev/null
echo "[3/4] detecting cSS candidates"
python $SC/detect_css.py $S -2.5 85 3
if [ -z "$NOVIEW" ]; then
  echo "[4/4] opening viewer"
  freeview -v $B/data/${S}_swi.nii \
    $B/review/${S}_candidates.nii.gz:colormap=lut:opacity=0.6 > /dev/null 2>&1 &
fi
echo ""
echo "NEXT: python $SC/export_review.py $S 40 | pbcopy    (paste into the workbook)"
echo "      review in freeview, then: python $SC/mark_css.py $S <cSS ids> [--ich <ids>]"
