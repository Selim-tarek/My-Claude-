#!/bin/zsh
# Bring T1 (and FLAIR) anatomy onto the SWI grid (v4).
# usage: prep_anat.sh SUBJECT_ID --t1 T1.nii(.gz) [--flair FLAIR.nii(.gz)] [--recon]
#   --recon  full FreeSurfer recon-all on the T1 (about 1-3 h on an M3): gives Destrieux sulcal
#            labels -> true SULCAL multifocality scoring. Re-used if subjects/ID already finished.
#   default  SynthSeg --parc --robust on the T1 (minutes): sharper cortex/CSF than SynthSeg on SWI,
#            but no sulcal labels (scoring stays Euclidean).
# Needs data/ID_swi.nii (run_css.sh imports it). Writes, all on the SWI grid, in work/:
#   ID_t1seg_swispace.nii.gz    aseg + DK cortex labels  -> align_seg.py / detect_css.py
#   ID_a2009s_swispace.nii.gz   Destrieux labels         -> score_css.py (recon-all only)
#   ID_flair_swispace.nii.gz    FLAIR                    -> detect_css.py (acute cSAH check)
#   ID_swi2t1.lta, ID_flair2swi.lta   registrations: ALWAYS check them with the printed freeview line
# Use the 3-D T1 (MPRAGE / SPGR / BRAVO) and the 3-D or 2-D FLAIR. Do NOT use the minIP series.
set -e
S=$1; shift 2>/dev/null || true
T1=""; FL=""; RECON=0
while [ $# -gt 0 ]; do
  case "$1" in
    --t1) T1=$2; shift 2;;
    --flair) FL=$2; shift 2;;
    --recon) RECON=1; shift;;
    *) echo "prep_anat.sh: unknown option $1"; exit 1;;
  esac
done
if [ -z "$S" ] || [ -z "$T1" ]; then
  echo "usage: prep_anat.sh SUBJECT_ID --t1 T1.nii [--flair FLAIR.nii] [--recon]"; exit 1
fi
B=${CSS_BASE:-$HOME/css_project}; W=$B/work
# always the project folder (FreeSurfer's setup script points SUBJECTS_DIR at its own install dir);
# set CSS_SUBJECTS_DIR to use a different folder
export SUBJECTS_DIR=${CSS_SUBJECTS_DIR:-$B/subjects}
SWI=$B/data/${S}_swi.nii
[ -f $SWI ] || { echo "ERROR: $SWI not found - run run_css.sh $S <swi> first"; exit 1; }
mkdir -p $W $SUBJECTS_DIR

echo "[a] T1: import + skull-strip (SynthStrip)"
mri_convert "$T1" $B/data/${S}_t1.nii.gz > /dev/null
mri_synthstrip -i $B/data/${S}_t1.nii.gz -o $W/${S}_t1_brain.nii.gz > /dev/null

A2=""
if [ $RECON -eq 1 ] || [ -f $SUBJECTS_DIR/$S/mri/aparc.a2009s+aseg.mgz ]; then
  if [ ! -f $SUBJECTS_DIR/$S/mri/aparc.a2009s+aseg.mgz ]; then
    echo "[b] recon-all -all (1-3 h; log: $W/${S}_recon.log)"
    recon-all -s $S -i $B/data/${S}_t1.nii.gz -all -threads 8 > $W/${S}_recon.log 2>&1
  else
    echo "[b] recon-all already finished for $S - re-using it"
  fi
  LAB=$SUBJECTS_DIR/$S/mri/aparc+aseg.mgz
  A2=$SUBJECTS_DIR/$S/mri/aparc.a2009s+aseg.mgz
  echo "[c] register SWI -> T1 (bbregister, boundary-based, T2*-like contrast)"
  bbregister --s $S --mov $SWI --reg $W/${S}_swi2t1.lta --t2 --init-coreg > $W/${S}_bbreg.log 2>&1
  grep -i "min cost" $W/${S}_bbreg.log | tail -1 || true
  REGM=bbregister; COST=$(awk 'NR==1{print $1}' $W/${S}_swi2t1.lta.mincost 2>/dev/null || true)
else
  echo "[b] SynthSeg --parc --robust on the T1"
  mri_synthseg --i $B/data/${S}_t1.nii.gz --o $W/${S}_t1synthseg.nii.gz --parc --robust --threads 8 > /dev/null
  LAB=$W/${S}_t1synthseg.nii.gz
  echo "[c] register SWI -> T1 (mri_coreg, rigid, mutual information)"
  mri_coreg --mov $SWI --ref $W/${S}_t1_brain.nii.gz --reg $W/${S}_swi2t1.lta --dof 6 > $W/${S}_coreg.log 2>&1
  REGM=mri_coreg; COST=""
fi

# v4.14 QC: registration method + cost (bbregister mincost; mri_coreg reports none) -> review/ID_qc.json
echo "$COST" | grep -Eq '^[0-9.eE+-]+$' || COST=null
printf '{"method": "%s", "cost": %s}\n' "$REGM" "$COST" > $W/${S}_reg.json
echo "[d] labels -> SWI grid (nearest neighbour)"
mri_vol2vol --mov $SWI --targ $LAB --lta $W/${S}_swi2t1.lta --inv --interp nearest \
  --o $W/${S}_t1seg_swispace.nii.gz > /dev/null
if [ -n "$A2" ]; then
  mri_vol2vol --mov $SWI --targ $A2 --lta $W/${S}_swi2t1.lta --inv --interp nearest \
    --o $W/${S}_a2009s_swispace.nii.gz > /dev/null
fi

if [ -n "$FL" ]; then
  echo "[e] FLAIR: import, skull-strip, register -> SWI, resample"
  mri_convert "$FL" $B/data/${S}_flair.nii.gz > /dev/null
  mri_synthstrip -i $B/data/${S}_flair.nii.gz -o $W/${S}_flair_brain.nii.gz > /dev/null
  mri_coreg --mov $W/${S}_flair_brain.nii.gz --ref $SWI --reg $W/${S}_flair2swi.lta --dof 6 > $W/${S}_flaircoreg.log 2>&1
  mri_vol2vol --mov $W/${S}_flair_brain.nii.gz --targ $SWI --lta $W/${S}_flair2swi.lta \
    --o $W/${S}_flair_swispace.nii.gz > /dev/null
fi

echo ""
echo "QC - labels must hug the SWI cortex, FLAIR must line up (toggle layers):"
echo "  freeview -v $SWI $W/${S}_t1seg_swispace.nii.gz:colormap=lut:opacity=0.3" \
     "$( [ -n "$FL" ] && echo "$W/${S}_flair_swispace.nii.gz:visible=0" )" \
     "$( [ -n "$A2" ] && echo "$W/${S}_a2009s_swispace.nii.gz:colormap=lut:opacity=0.3:visible=0" )"
echo "If misaligned: delete work/${S}_t1seg_swispace.nii.gz and the pipeline falls back to SWI SynthSeg."
