#!/bin/bash
# Self-test with a synthetic phantom brain - no patient data, no FreeSurfer needed.
# usage: bash tests/run_tests.sh      (from the package root, inside the 'css' conda env)
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
T=$(mktemp -d); export CSS_BASE=$T FREESURFER_HOME=$T/fs CSS_REVIEW_TEST=1
mkdir -p $T/data $T/synthseg $T/work $T/review $T/results $T/fs $T/scripts
cp $ROOT/scripts/*.py $T/scripts/
printf '0 Unknown 0 0 0 0\n2 Left-Cerebral-White-Matter 245 245 245 0\n24 CSF 60 60 60 0\n41 Right-Cerebral-White-Matter 245 245 245 0\n1022 ctx-lh-postcentral 220 20 20 0\n1028 ctx-lh-superiorfrontal 20 220 160 0\n2022 ctx-rh-postcentral 220 20 20 0\n2028 ctx-rh-superiorfrontal 20 220 160 0\n' > $T/fs/FreeSurferColorLUT.txt
S=$T/scripts
echo "== phantoms";            python $ROOT/tests/phantom.py PH1 1; python $ROOT/tests/phantom.py PH2 2
echo "== align";               python $S/align_seg.py PH1; python $S/align_seg.py PH2
echo "== detect (clean)";      python $S/detect_css.py PH1 | head -1
echo "== export";              python $S/export_review.py PH1 3
echo "== synthetic lesions";   python $S/make_synthetic.py PH1 8 0.6 1 | head -1
cp $T/synthseg/PH1_seg.nii.gz $T/synthseg/PH1S_seg.nii.gz; python $S/align_seg.py PH1S > /dev/null
python $S/detect_css.py PH1S | head -1
R=$(python $S/eval_synthetic.py PH1S | grep RESULT); echo "$R"
FOUND=$(echo "$R" | sed -E 's/.*found=([0-9]+).*/\1/')
[ "$FOUND" -ge 4 ] || { echo "FAIL: sensitivity too low ($FOUND/8)"; exit 1; }
echo "== scoring with truth";  python $S/score_css.py PH1S --truth > $T/score_test.log; grep -A2 "cSS multifocality" $T/score_test.log
echo "== interactive reviewer (headless test mode)"; python $S/review_css.py PH1S --top 8 --redo > $T/review_test.log; grep "cSS multifocality" $T/review_test.log
echo "== feature report";      python $S/feature_report.py | head -3
echo "== stress test (small)"; python $S/stress_test.py --hosts PH1,PH2 --depths 0.6 --seeds 1 --tag test | grep -A2 SUMMARY
echo; echo "ALL TESTS PASSED  (temp dir $T)"
