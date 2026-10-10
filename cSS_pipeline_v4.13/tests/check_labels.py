"""v4.15 checks: expert sheet with excluded objects, fine-grained label import, long label table,
missed-lesion fates, study metadata with a fixed split.   usage: check_labels.py   (needs PH3 in $CSS_BASE)"""
import os, sys, subprocess, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
B = os.environ["CSS_BASE"]; S = f"{B}/scripts"; py = sys.executable
def fail(m): print("FAIL:", m); sys.exit(1)
def run(*a, ok=True):
    r = subprocess.run([py, *a], capture_output=True, text=True)
    if ok and r.returncode: print(r.stdout[-800:], r.stderr[-1500:]); fail("step failed: " + " ".join(a))
    return r
run(f"{S}/detect_css.py", "PH3")
truth = nib.load(f"{B}/work/PH3_truth.nii.gz").get_fdata() > 0
cand = nib.load(f"{B}/review/PH3_candidates.nii.gz").get_fdata().astype(int)
excl = nib.load(f"{B}/review/PH3_excluded.nii.gz").get_fdata().astype(int)
drop = nib.load(f"{B}/work/PH3_dropped.nii.gz").get_fdata().astype(int)
nx = len(pd.read_csv(f"{B}/review/PH3_excluded.csv"))

# 1. sheet mixes the excluded objects in, unmarked, and the key says which list each number is
print("   " + run(f"{S}/expert_sheet.py", "PH3", "--include-excluded").stdout.splitlines()[0])
key = pd.read_csv(f"{B}/review/PH3_expert_key.csv")
if (key.list == "excluded").sum() != nx or (key.list == "candidate").sum() != len(pd.read_csv(f"{B}/review/PH3_candidates.csv")):
    fail(f"key lists {key.list.value_counts().to_dict()}")
if not key.sheet_no.tolist() == list(range(1, len(key) + 1)): fail("sheet numbers")

# 2. fine codes + confidence; one number left out (-> N, defaulted); excluded veins labelled VC
ans, left_out = [], int(key.sheet_no.iloc[-1])
for _, r in key.iterrows():
    if r.sheet_no == left_out: continue
    m = (cand == r.obj_id) if r.list == "candidate" else (excl == r.obj_id)
    ans.append(f"{r.sheet_no}{'C:5' if truth[m].mean() > 0.2 else ('VS:3' if r.list == 'candidate' else 'VC:4')}")
out = run(f"{S}/expert_import.py", "PH3", "expertB", "--letters", " ".join(ans)).stdout
print("   " + out.strip().replace("\n", "\n   "))
e = pd.read_csv(f"{B}/review/PH3_expert_expertB.csv")
if list(e.columns[:6]) != ["cand_id", "call", "vox_i", "vox_j", "vox_k", "volume_mm3"]: fail("review expert file columns")
if not set(e.call) <= {"cSS", "Vein", "Normal"} or not {"cSS"} <= set(e.call): fail(f"coarse calls {set(e.call)}")
L = pd.read_csv(f"{B}/labels/candidate_labels.csv"); Lb = L[L.reader_id == "expertB"]
if (Lb.list == "excluded").sum() != nx: fail("excluded objects missing from the label table")
if Lb.defaulted.sum() != 1 or Lb[Lb.sheet_no == left_out].label.iloc[0] != "normal_other": fail("default N")
if not set(Lb[Lb.list == "excluded"].label) <= {"vein_sulcal", "cSS"}: fail("excluded labels")
if Lb[Lb.label == "cSS"].confidence.ne(5).any(): fail("confidence not stored")
# re-import of the same session replaces, a second session adds
run(f"{S}/expert_import.py", "PH3", "expertB", "--letters", " ".join(ans))
if len(pd.read_csv(f"{B}/labels/candidate_labels.csv").query("reader_id == 'expertB'")) != len(Lb): fail("re-import duplicated rows")
xn = int(key[key.list == "excluded"].sheet_no.iloc[0])
out2 = run(f"{S}/expert_import.py", "PH3", "expertB", "--session", "2", "--letters",
           " ".join(a for a in ans if not a.startswith(f"{xn}")) + f" {xn}C").stdout
if "EXCLUDED object(s) cSS" not in out2: fail("no warning when an excluded object is called cSS")
L2 = pd.read_csv(f"{B}/labels/candidate_labels.csv").query("reader_id == 'expertB'")
if sorted(L2.session.unique()) != [1, 2] or len(L2) != 2 * len(Lb): fail("session 2 rows")
if run(f"{S}/expert_import.py", "PH3", "expertB", "--letters", "1QQ", ok=False).returncode == 0: fail("bad code accepted")
# old key format (no list/obj_id columns) still imports
key[key.list == "candidate"].drop(columns=["list", "obj_id"]).to_csv(f"{B}/review/PH3_expert_key.csv", index=False)
run(f"{S}/expert_import.py", "PH3", "expertC", "--css", ",".join(str(n) for n, a in
    zip(key.sheet_no, ans) if a.endswith("C:5") and n in set(key[key.list == "candidate"].sheet_no)))
print("   expert codes, confidence, defaults, sessions, excluded objects and the old key format: OK")

# 3. missed-lesion outlines: one blob per fate
vox = np.array(nib.load(f"{B}/data/PH3_swi.nii").header.get_zooms()[:3])
M = np.zeros(cand.shape, np.uint8)
def piece(obj):
    """one connected outline: the object's voxels in a 3x3x3 box around one of its voxels"""
    q = np.argwhere(obj); q = q[len(q) // 2]; box = np.zeros(obj.shape, bool)
    box[q[0] - 1:q[0] + 2, q[1] - 1:q[1] + 2, q[2] - 1:q[2] + 2] = True
    return obj & box
far_from = ndi.distance_transform_edt((cand == 0) & (excl == 0), sampling=vox)
M[piece((cand == 1) & truth)] = 1                                   # shown
M[piece((excl > 0) & (ndi.distance_transform_edt(cand == 0, sampling=vox) > 4))] = 1    # excluded
d_ok = [i for i in np.unique(drop) if i and far_from[drop == i].min() > 4]
M[piece(drop == d_ok[0])] = 1                                       # dropped (gate reason)
seg = nib.load(f"{B}/work/PH3_seg_swispace.nii.gz").get_fdata().astype(int)
free = np.isin(seg, [2, 41]) & (ndi.distance_transform_edt((cand == 0) & (excl == 0) & (drop == 0), sampling=vox) > 6)
p = np.argwhere(free)[len(np.argwhere(free)) // 2]
M[p[0] - 1:p[0] + 2, p[1] - 1:p[1] + 2, p[2]:p[2] + 2] = 1          # nothing there -> not generated
os.makedirs(f"{B}/labels", exist_ok=True)
aff = nib.load(f"{B}/data/PH3_swi.nii").affine
nib.save(nib.Nifti1Image(M, aff), f"{B}/labels/PH3_missed_expertB.nii.gz")
print("   " + run(f"{S}/label_missed.py", "PH3", "expertB").stdout.strip().replace("\n", "\n   "))
ml = pd.read_csv(f"{B}/labels/missed_lesions.csv")
if sorted(ml.fate) != sorted(["shown", "excluded", "dropped", "not generated"]): fail(f"fates {ml.fate.tolist()}")
if not ml[ml.fate == "excluded"].rule_version.iloc[0].startswith("E"): fail("excluded fate without rule id")
if not str(ml[ml.fate == "dropped"].reason.iloc[0]).startswith(("volume", "<2", "elong", "no cortex")): fail("dropped reason")
run(f"{S}/label_missed.py", "PH3", "expertB")
if len(pd.read_csv(f"{B}/labels/missed_lesions.csv")) != 4: fail("re-run duplicated missed-lesion rows")

# 4. study metadata: QC fields copied, split fixed once assigned
print("   " + run(f"{S}/study_meta.py", "PH3", "--status", "pos", "--ref-L", "1", "--ref-R", "2", "--split", "dev").stdout.strip())
if run(f"{S}/study_meta.py", "PH3", "--split", "test", ok=False).returncode == 0: fail("split changed without --force")
run(f"{S}/study_meta.py", "PH3", "--motion", "1")
sm = pd.read_csv(f"{B}/labels/study_meta.csv").set_index("study_id").loc["PH3"]
if sm.split != "dev" or sm.css_status != "pos" or int(sm.ref_multifocality_0_4) != 3 or sm.qc_status != "low_confidence_review":
    fail(f"study_meta row {sm.to_dict()}")
run(f"{S}/study_meta.py", "PH3", "--split", "test", "--force")
if "forced" not in pd.read_csv(f"{B}/labels/study_meta.csv").set_index("study_id").loc["PH3", "note"]: fail("forced change not logged")
print("   study_meta: QC fields copied, split fixed unless --force (logged)")
print("RESULT labels OK")
