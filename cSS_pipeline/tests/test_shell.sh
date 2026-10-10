#!/bin/bash
# Shell-flow test for run_css.sh + prep_anat.sh with STUBBED FreeSurfer tools and python
# (checks argument parsing, which tools run in which order, and the re-use / staleness logic).
# usage: bash tests/test_shell.sh        Uses real zsh if installed, else bash as a stand-in.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
T=$(mktemp -d); mkdir -p $T/bin $T/base/scripts; cp $ROOT/scripts/*.py $ROOT/scripts/*.sh $T/base/scripts/
cat > $T/bin/fsstub <<'EOS'
#!/bin/bash
n=$(basename $0); echo "$n $*" >> $STUBLOG
out=""; in=""; targ=""
while [ $# -gt 0 ]; do case "$1" in
  -o|--o) out=$2; shift 2;; -i|--i|--mov) [ -z "$in" ] && in=$2; shift 2;; --targ) targ=$2; shift 2;;
  --reg) touch $2; shift 2;;
  -s) [ "$n" = recon-all ] && mkdir -p $SUBJECTS_DIR/$2/mri && touch $SUBJECTS_DIR/$2/mri/aparc+aseg.mgz $SUBJECTS_DIR/$2/mri/aparc.a2009s+aseg.mgz; shift 2;;
  *) [ "$n" = mri_convert ] && { [ -z "$in" ] && in=$1 || out=$1; }; shift;; esac; done
[ -n "$out" ] && [ -n "$in" ] && cp "$in" "$out" 2>/dev/null || true
EOS
chmod +x $T/bin/fsstub
for t in mri_convert mri_synthstrip mri_synthseg mri_coreg mri_vol2vol bbregister recon-all freeview; do ln -s fsstub $T/bin/$t; done
printf '#!/bin/bash\necho "python $(basename $1) ${*:2}" >> $STUBLOG\n' > $T/bin/python; chmod +x $T/bin/python
command -v zsh >/dev/null || { printf '#!/bin/bash\nexec bash "$@"\n' > $T/bin/zsh; chmod +x $T/bin/zsh; echo "(zsh not installed: using bash as a stand-in)"; }
echo swi > $T/swi.nii; echo t1 > $T/t1.nii; echo fl > $T/fl.nii
export PATH=$T/bin:$PATH STUBLOG=$T/log CSS_BASE=$T/base NOVIEW=1 SUBJECTS_DIR=/somewhere/else  # must be ignored
SH="zsh"
fail() { echo "FAIL: $1"; cat $T/log; exit 1; }
seen() { grep -q "$1" $T/log || fail "expected call: $1"; }
notseen() { ! grep -q "$1" $T/log || fail "unexpected call: $1"; }

echo "== 1. SWI only"
: > $T/log; $SH $T/base/scripts/run_css.sh P1 $T/swi.nii > /dev/null
seen "mri_synthstrip -i $T/base/work/P1_raw"; seen "mri_synthseg --i $T/base/data/P1_swi.nii"
seen "python align_seg.py P1"; seen "python detect_css.py P1 -2.5 85 3"; notseen mri_coreg
echo "== 2. re-run: SynthSeg re-used"
: > $T/log; $SH $T/base/scripts/run_css.sh P1 > /dev/null; notseen mri_synthseg; seen "detect_css.py P1"
echo "== 3. SWI re-imported -> stale SynthSeg is redone"
sleep 1; : > $T/log; $SH $T/base/scripts/run_css.sh P1 $T/swi.nii > /dev/null; seen "mri_synthseg --i $T/base/data/P1_swi.nii"
echo "== 4. SWI + T1 + FLAIR (SynthSeg on T1, mri_coreg)"
: > $T/log; $SH $T/base/scripts/run_css.sh P2 $T/swi.nii --t1 $T/t1.nii --flair $T/fl.nii > /dev/null
seen "mri_synthseg --i $T/base/data/P2_t1.nii.gz"; seen "mri_coreg --mov $T/base/data/P2_swi.nii --ref $T/base/work/P2_t1_brain"
seen "mri_vol2vol --mov $T/base/data/P2_swi.nii --targ $T/base/work/P2_t1synthseg.nii.gz --lta $T/base/work/P2_swi2t1.lta --inv --interp nearest"
seen "mri_coreg --mov $T/base/work/P2_flair_brain.nii.gz --ref $T/base/data/P2_swi.nii"
seen "mri_vol2vol --mov $T/base/work/P2_flair_brain.nii.gz --targ $T/base/data/P2_swi.nii"
notseen "mri_synthseg --i $T/base/data/P2_swi.nii"; notseen recon-all; notseen bbregister
echo "== 5. --recon (recon-all, bbregister, Destrieux)"
: > $T/log; $SH $T/base/scripts/run_css.sh P3 $T/swi.nii --t1 $T/t1.nii --recon > /dev/null
seen "recon-all -s P3"; seen "bbregister --s P3 --mov $T/base/data/P3_swi.nii --reg $T/base/work/P3_swi2t1.lta --t2 --init-coreg"
seen "aparc.a2009s+aseg.mgz --lta $T/base/work/P3_swi2t1.lta --inv --interp nearest --o $T/base/work/P3_a2009s_swispace"
echo "== 6. recon-all already done -> re-used without --recon"
: > $T/log; $SH $T/base/scripts/prep_anat.sh P3 --t1 $T/t1.nii > /dev/null; notseen "recon-all"; seen bbregister
echo "== 6b. --phase is imported for the reviewer"
echo ph > $T/ph.nii; : > $T/log; $SH $T/base/scripts/run_css.sh P4 $T/swi.nii --phase $T/ph.nii > /dev/null
seen "mri_convert $T/ph.nii $T/base/data/P4_phase.nii.gz"
echo "== 7. batch run_all.sh reaches detection for every subject"
: > $T/log; $SH $T/base/scripts/run_all.sh P1 P2 > /dev/null || true   # stub python prints no summary
[ $(grep -c "detect_css.py" $T/log) -eq 2 ] || fail "run_all.sh did not run detection for both subjects"
echo "SHELL TESTS PASSED"
