"""v4.17 checks: Cohen's / weighted kappa, two-reader consensus with adjudication (consensus.py).
usage: check_consensus.py    (needs PH3 in $CSS_BASE)"""
import os, sys, json, subprocess, numpy as np, nibabel as nib, pandas as pd
B = os.environ["CSS_BASE"]; S = f"{B}/scripts"; py = sys.executable
sys.path.insert(0, S)
from css_common import cohen_kappa
def fail(m): print("FAIL:", m); sys.exit(1)
def run(*a):
    r = subprocess.run([py, *a], capture_output=True, text=True)
    if r.returncode: print(r.stdout[-800:], r.stderr[-1500:]); fail("step failed: " + " ".join(a))
    return r.stdout

# 1. kappa: textbook values
if abs(cohen_kappa([0, 0, 1, 1], [0, 1, 1, 1]) - 0.5) > 1e-9: fail("Cohen's kappa")
if abs(cohen_kappa([0, 1, 2, 3, 4], [0, 1, 2, 3, 4], range(5), "quadratic") - 1) > 1e-9: fail("weighted kappa = 1")
w1 = cohen_kappa([0, 1, 2, 3, 4], [1, 2, 3, 4, 4], range(5), "quadratic")
w2 = cohen_kappa([0, 1, 2, 3, 4], [4, 3, 2, 1, 0], range(5), "quadratic")
if not (0.75 < w1 < 1 and w2 < -0.9): fail(f"weighted kappa ordering {w1:.2f} {w2:.2f}")
print(f"   kappa: textbook 0.50, quadratic-weighted near-miss {w1:.2f}, reversed {w2:.2f}")

# 2. two readers on PH3 (sheet with excluded objects): R2 disagrees on two classes and one fine code
run(f"{S}/detect_css.py", "PH3"); run(f"{S}/expert_sheet.py", "PH3", "--include-excluded")
key = pd.read_csv(f"{B}/review/PH3_expert_key.csv")
tr = nib.load(f"{B}/work/PH3_truth.nii.gz").get_fdata() > 0
C = nib.load(f"{B}/review/PH3_candidates.nii.gz").get_fdata().astype(int)
X = nib.load(f"{B}/review/PH3_excluded.nii.gz").get_fdata().astype(int)
r1 = {}
for _, r in key.iterrows():
    m = (C == r.obj_id) if r.list == "candidate" else (X == r.obj_id)
    r1[int(r.sheet_no)] = "C" if tr[m].mean() > 0.2 else "VC"
nums = sorted(r1); c_ = [n for n in nums if r1[n] == "C"]; v_ = [n for n in nums if r1[n] == "VC"]
r2 = dict(r1); r2[c_[0]] = "VS"; r2[v_[0]] = "A"; r2[v_[1]] = "VT"          # 2 class disagreements, 1 fine-code
for rd, ans in (("readerR1", r1), ("readerR2", r2)):
    run(f"{S}/expert_import.py", "PH3", rd, "--letters", " ".join(f"{n}{c}" for n, c in ans.items()))
out = run(f"{S}/consensus.py", "--readers", "readerR1,readerR2", "--subjects", "PH3")
print("   " + out.strip().replace("\n", "\n   "))
todo = pd.read_csv(f"{B}/labels/adjudication_todo.csv")
if len(todo) != 2: fail(f"{len(todo)} pending, expected the 2 class disagreements")
cons = pd.read_csv(f"{B}/labels/consensus.csv"); cons = cons[cons.study_id == "PH3"]
fine = cons[(cons.method == "agreement") & (cons.consensus_label == "other_vein")]
if len(fine) != 1: fail("VC vs VT should agree at class level (other_vein)")
st = json.load(open(f"{B}/labels/agreement_stats.json"))
if not (0 < st["estimates"]["kappa_class"] < 1): fail(f"kappa {st['estimates']}")
# 3. adjudication: a third reader resolves both, --apply closes them
todo["adjudicated_label"] = ["C", "VC"]; todo["adjudicator"] = "readerR3"
todo.to_csv(f"{B}/labels/adjudication_todo.csv", index=False)
run(f"{S}/consensus.py", "--readers", "readerR1,readerR2", "--apply", f"{B}/labels/adjudication_todo.csv")
cons = pd.read_csv(f"{B}/labels/consensus.csv"); cons = cons[cons.study_id == "PH3"]
if (cons.method == "pending").any() or (cons.method == "adjudicated").sum() != 2: fail("adjudication not applied")
if len(cons) != len(key): fail("consensus does not cover every object on the sheet")
print(f"   consensus: {cons.method.value_counts().to_dict()}; kappa(6 classes) {st['estimates']['kappa_class']}")
print("RESULT consensus OK")
