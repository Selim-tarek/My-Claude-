#!/bin/bash
# Self-test with a synthetic phantom brain - no patient data, no FreeSurfer needed.
# usage: bash tests/run_tests.sh      (from the package root, inside the 'css' conda env)
# Shell flow of run_css.sh / prep_anat.sh with stubbed FreeSurfer: bash tests/test_shell.sh
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
echo "== calls survive a detector re-run (fingerprint matching)"
cp $T/review/PH1S_candidates.csv $T/old_cand.csv
python $S/detect_css.py PH1S -2.3 80 3 > /dev/null
python - "$T" <<'PYEOF'
import sys, pandas as pd; sys.path.insert(0, sys.argv[1] + "/scripts")
from css_common import load_calls, FP
b = sys.argv[1]; old = pd.read_csv(b + "/old_cand.csv"); new = pd.read_csv(b + "/review/PH1S_candidates.csv")
raw = pd.read_csv(b + "/review/PH1S_calls.csv"); calls, warn = load_calls(b, "PH1S", new)
for cid, call in calls.items():          # every re-attached call must sit on the same lesion as before
    r = raw[raw.call == call]; n = new.set_index("cand_id").loc[cid]
    assert (((r[FP[:3]] - n[FP[:3]].astype(float)) ** 2).sum(axis=1) ** .5 <= 2).any(), (cid, call)
print(f"   {len(calls)}/{len(raw)} calls re-attached by fingerprint" + (f"; {warn}" if warn else ""))
PYEOF
python $S/detect_css.py PH1S > /dev/null
echo "== feature report";      python $S/feature_report.py | head -3
echo "== stress test (small)"; python $S/stress_test.py --hosts PH1,PH2 --depths 0.6 --seeds 1 --tag test | grep -A2 SUMMARY
echo "== v4 phantom: tram-track / convexity cSS vs tubular & surface veins, ICH, sulcal scoring"
python $ROOT/tests/phantom_v4.py PH3 1; python $S/align_seg.py PH3 > /dev/null
python $S/detect_css.py PH3 | head -1
python $ROOT/tests/check_v4.py PH3 | tail -3
python $S/score_css.py PH3 --truth > $T/score_v4.log
grep -q "3/4" $T/score_v4.log && grep -q "sulcal (Destrieux)" $T/score_v4.log || { cat $T/score_v4.log; echo "FAIL: sulcal score on PH3 should be 3/4 (L1 + R2)"; exit 1; }
grep -A3 "cSS multifocality" $T/score_v4.log
echo "== reviewer with a phase image (calcium vs blood-product panel)"
python - "$T" <<'PYEOF'
import sys, numpy as np, nibabel as nib
b = sys.argv[1]; r = nib.load(b + "/data/PH3_swi.nii")
ph = np.random.default_rng(0).integers(-4096, 4096, r.shape).astype(np.int16)   # scanner-unit phase
nib.save(nib.Nifti1Image(ph, r.affine), b + "/data/PH3_phase.nii.gz")
PYEOF
CSS_REVIEW_TEST=1 python $S/review_css.py PH3 --top 2 --redo > /dev/null && [ -s $T/review/PH3_review_preview.png ] && echo "   phase panel rendered"
rm $T/data/PH3_phase.nii.gz
python $S/score_css.py PH3 --truth > $T/score_v4.log
echo "== ICH rule in sulci (Charidimou 2017): reader-drawn ICH -> left foci next to it excluded"
cp $T/work/PH3_ich_truth.nii.gz $T/work/PH3_ich.nii.gz
python $S/score_css.py PH3 --truth > $T/score_ich.log; rm $T/work/PH3_ich.nii.gz
grep -q "2/4" $T/score_ich.log && grep -q "EXCLUDED" $T/score_ich.log || { cat $T/score_ich.log; echo "FAIL: ICH sulcal rule"; exit 1; }
grep "ICH rule" $T/score_ich.log
echo "== Boston v2.0 cSS count by gyri: all foci -> >=2; left only (2 adjacent gyri) -> 1"
grep -q "Boston v2.0 cSS count: >=2" $T/score_v4.log || { cat $T/score_v4.log; echo "FAIL: Boston v2.0 >=2"; exit 1; }
python -c "import pandas as pd,sys; p=sys.argv[1]+'/review/PH3_candidates.csv'; d=pd.read_csv(p); d.loc[d.hemi=='R','accept']=0; d.to_csv(p,index=False)" $T
python $S/score_css.py PH3 > $T/score_b2.log
grep -q "Boston v2.0 cSS count: 1 focus" $T/score_b2.log || { cat $T/score_b2.log; echo "FAIL: Boston v2.0 single focus"; exit 1; }
grep "Boston" $T/score_v4.log $T/score_b2.log
echo "== expert sheet (blinded PDF) + import of the expert's letters"
python $S/detect_css.py PH3 > /dev/null
python $S/expert_sheet.py PH3 | head -1
[ -s $T/review/PH3_expert_sheet.pdf ] || { echo "FAIL: no expert sheet"; exit 1; }
python - "$T" <<'PYEOF'
import sys, pandas as pd, nibabel as nib, numpy as np
b = sys.argv[1]; key = pd.read_csv(b + "/review/PH3_expert_key.csv")
truth = nib.load(b + "/work/PH3_truth.nii.gz").get_fdata() > 0
cand = nib.load(b + "/review/PH3_candidates.nii.gz").get_fdata().astype(int)
letters = " ".join(f"{n}{'C' if (truth[cand == c]).mean() > 0.2 else 'V'}" for n, c in zip(key.sheet_no, key.cand_id))
open(b + "/letters.txt", "w").write(letters)
PYEOF
python $S/expert_import.py PH3 expertA --letters "$(cat $T/letters.txt)" --score | grep -E "expert|multifocality"
python - "$T" <<'PYEOF'
import sys, pandas as pd, nibabel as nib
b = sys.argv[1]; e = pd.read_csv(b + "/review/PH3_expert_expertA.csv")
want = open(b + "/letters.txt").read().count("C")
assert (e.call == "cSS").sum() == want, (e.call.value_counts(), want)   # every C mapped back to its candidate
print("   expert calls mapped back to the right candidates")
PYEOF
python $S/feature_report.py | head -2
echo "== normal anatomy: a dark line on the cerebellum must be excluded, cSS lesions kept"
python $ROOT/tests/phantom_infra.py; python $S/align_seg.py PH4 > /dev/null
python $S/detect_css.py PH4 | head -1
python $S/check_known.py PH4 --old $T/work/PH3_truth.nii.gz --ids 1,2,3,4,5,6 | tail -1 | grep -q "6 of 6" \
  || { python $S/check_known.py PH4 --old $T/work/PH3_truth.nii.gz --ids 1,2,3,4,5,6; echo "FAIL: cSS lost on PH4"; exit 1; }
python $S/check_known.py PH4 --old $T/work/PH4_infra_line.nii.gz --ids 1 | head -1
python $S/check_known.py PH4 --old $T/work/PH4_infra_line.nii.gz --ids 1 | head -1 | grep -q "kept" \
  && { echo "FAIL: cerebellar line still shown as a candidate"; exit 1; }
echo "== partial-coverage slab (bottom cut off) with T1 anatomy: skull-base rules must switch off"
python - "$T" <<'PYEOF'
import sys, numpy as np, nibabel as nib
b = sys.argv[1]; img = nib.load(b + "/data/PH3_swi.nii"); I = img.get_fdata().astype(np.float32)
seg = nib.load(b + "/work/PH3_seg_swispace.nii.gz").get_fdata().astype(np.int32)
I[:, :, :30] = 0; seg[:, :, :30] = 0                      # slab starts in the middle of the brain
nib.save(nib.Nifti1Image(I, img.affine), b + "/data/PH5_swi.nii")
nib.save(nib.Nifti1Image(seg, img.affine), b + "/synthseg/PH5_seg.nii.gz")
nib.save(nib.Nifti1Image(seg, img.affine), b + "/work/PH5_t1seg_swispace.nii.gz")   # pretend T1 labels
PYEOF
python $S/align_seg.py PH5 > /dev/null
python $S/detect_css.py PH5 | head -2 | tee $T/ph5.log
grep -q "anatomy T1" $T/ph5.log && grep -q "cut off at the bottom" $T/ph5.log || { echo "FAIL: coverage / anatomy detection"; exit 1; }
grep -q "skull base" $T/ph5.log && { echo "FAIL: skull-base rule active on a cut-off slab"; exit 1; }
echo "== T1 anatomy on the cerebellum phantom: position rules active, cSS kept"
cp $T/synthseg/PH4_seg.nii.gz $T/work/PH4_t1seg_swispace.nii.gz; python $S/align_seg.py PH4 > /dev/null
python $S/detect_css.py PH4 | head -2
python $S/check_known.py PH4 --old $T/work/PH3_truth.nii.gz --ids 1,2,3,4,5,6 | tail -1 | grep -q "6 of 6" || { echo "FAIL: cSS lost (T1 mode)"; exit 1; }
python $S/check_known.py PH4 --old $T/work/PH4_infra_line.nii.gz --ids 1 | head -1 | grep -q "kept" && { echo "FAIL: cerebellar line kept (T1 mode)"; exit 1; }
echo "== stress test on a host with T1 anatomy (PH4)"
python $S/stress_test.py --hosts PH4 --depths 0.6 --seeds 1 --tag t1host | grep -E "anatomy|darkness"
echo; echo "ALL TESTS PASSED  (temp dir $T)"
