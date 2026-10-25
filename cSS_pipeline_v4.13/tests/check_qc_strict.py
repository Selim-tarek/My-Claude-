"""v5.0 checks: CSS_QC_STRICT=1 (detect_css.py, align_seg.py, rank_candidates.py). Default (unset) = v4.19 behaviour.
usage: check_qc_strict.py    (needs PH3 from phantom_v4.py, aligned, in $CSS_BASE)
1. insufficient_quality (labels cover too few brain slices): strict -> header-only candidates CSV, exit 3;
   unset -> candidates written as before.
2. low_confidence_review (SWI-only anatomy): strict -> detection runs (exit 0) and rank_candidates.py multiplies
   score_model by 0.8; unset -> score_model unchanged.
3. align_seg.py: T1 labels older than the SWI -> strict stops with an error; unset -> silent SynthSeg fall-back."""
import os, sys, json, time, subprocess, numpy as np, nibabel as nib, pandas as pd
B = os.environ["CSS_BASE"]; S = f"{B}/scripts"; py = sys.executable
def fail(m): print("FAIL:", m); sys.exit(1)
def run(*a, strict=False, env=None):
    e = dict(os.environ); e.pop("CSS_QC_STRICT", None); e.update(env or {})
    if strict: e["CSS_QC_STRICT"] = "1"
    return subprocess.run([py, *a], capture_output=True, text=True, env=e)
img = nib.load(f"{B}/data/PH3_swi.nii"); I = img.get_fdata().astype(np.float32); aff = img.affine
seg = nib.load(f"{B}/synthseg/PH3_seg.nii.gz").get_fdata().astype(np.int32)
def case(name, sg):
    nib.save(nib.Nifti1Image(I, aff), f"{B}/data/{name}_swi.nii")
    nib.save(nib.Nifti1Image(sg, aff), f"{B}/synthseg/{name}_seg.nii.gz")
    p = f"{B}/work/{name}_t1seg_swispace.nii.gz"
    if os.path.exists(p): os.remove(p)
    r = run(f"{S}/align_seg.py", name)
    if r.returncode: fail(f"align_seg {name}: {r.stderr[-300:]}")

# 1. labels on the lower 40 % of slices only -> label_slice_coverage < 0.9 -> insufficient_quality
sg = seg.copy(); sg[:, :, int(0.4 * sg.shape[2]):] = 0
case("PHQ1", sg)
r0 = run(f"{S}/detect_css.py", "PHQ1")
q = json.load(open(f"{B}/review/PHQ1_qc.json"))
if r0.returncode or q["status"] != "insufficient_quality": fail(f"unset: exit {r0.returncode}, QC {q['status']}")
n0 = len(pd.read_csv(f"{B}/review/PHQ1_candidates.csv"))
r1 = run(f"{S}/detect_css.py", "PHQ1", strict=True)
c1 = pd.read_csv(f"{B}/review/PHQ1_candidates.csv")
if r1.returncode != 3 or len(c1) or "cannot be scored" not in r1.stdout: fail(f"strict: exit {r1.returncode}, {len(c1)} rows")
if "region_z" not in c1.columns: fail("strict empty CSV lost the column contract")
print(f"insufficient_quality ({q['reasons'][0]}): unset -> {n0} candidates, exit 0; strict -> 0 rows, exit 3, columns kept")

# 2. low confidence (SWI-only anatomy) + a small hand-made model
case("PHQ2", seg)
r = run(f"{S}/detect_css.py", "PHQ2", strict=True)
q = json.load(open(f"{B}/review/PHQ2_qc.json"))
if r.returncode or q["status"] != "low_confidence_review": fail(f"low-confidence case: exit {r.returncode}, {q['status']}")
model = dict(features=["darkness_z", "tube_ratio"], weights=[1.0, -1.0], bias=-2.0, impute_median=[3.0, 0.5],
             mean=[3.0, 0.5], sd=[1.0, 0.3], training_studies=["XX"], n_patients=3, lopo_auc=None, baseline="score_v4")
json.dump(model, open(f"{B}/work/qc_model.json", "w"))
def ranked(strict):
    rr = run(f"{S}/rank_candidates.py", "PHQ2", "--model", f"{B}/work/qc_model.json", strict=strict)
    if rr.returncode: fail(f"rank_candidates: {rr.stderr[-300:]}")
    return pd.read_csv(f"{B}/review/PHQ2_candidates.csv")
a, b = ranked(False), ranked(True)
ratio = (b.score_model / a.score_model).round(3).unique()
if len(ratio) != 1 or abs(ratio[0] - 0.8) > 0.002 or not (a.rank_model == b.rank_model).all():
    fail(f"strict low-confidence penalty: ratios {ratio}")
print(f"low_confidence_review: strict detection runs (exit 0); score_model x{ratio[0]} with the order unchanged")

# 3. stale T1 labels
nib.save(nib.Nifti1Image(seg, aff), f"{B}/work/PHQ2_t1seg_swispace.nii.gz")
old = time.time() - 3600; os.utime(f"{B}/work/PHQ2_t1seg_swispace.nii.gz", (old, old))
r0, r1 = run(f"{S}/align_seg.py", "PHQ2"), run(f"{S}/align_seg.py", "PHQ2", strict=True)
if r0.returncode or "synthseg" not in r0.stdout.lower() and "PHQ2_seg" not in r0.stdout: fail(f"unset: {r0.stdout} {r0.stderr}")
if r1.returncode == 0 or "older than the SWI" not in r1.stderr: fail(f"strict stale T1 labels accepted: {r1.stdout}")
print("stale T1 labels: unset -> SynthSeg fall-back as before; strict -> stops with an error")
print("RESULT qc_strict OK")
