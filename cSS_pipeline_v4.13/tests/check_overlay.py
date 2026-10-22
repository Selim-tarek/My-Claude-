"""v4.16 checks for surface_overlay.py on PH6 (needs phantom_surfaces.py PH6 and check_surfvein.py first).
Colour by detector list and by expert calls; annotation / id overlay consistency; Destrieux sulcus per object;
a lower-priority class never hides a cSS call."""
import os, sys, subprocess, numpy as np, nibabel as nib, pandas as pd
from nibabel.freesurfer.io import read_annot
B = os.environ["CSS_BASE"]; S = f"{B}/scripts"; py = sys.executable
def fail(m): print("FAIL:", m); sys.exit(1)
def run(*a, ok=True):
    r = subprocess.run([py, *a], capture_output=True, text=True)
    if ok and r.returncode: print(r.stdout[-800:], r.stderr[-1500:]); fail("step failed: " + " ".join(a))
    return r
tr = nib.load(f"{B}/work/PH6_truth.nii.gz").get_fdata().astype(int)
sv = nib.load(f"{B}/work/PH6_surfvein.nii.gz").get_fdata().astype(int)
C = nib.load(f"{B}/review/PH6_candidates.nii.gz").get_fdata().astype(int)
X = nib.load(f"{B}/review/PH6_excluded.nii.gz").get_fdata().astype(int)

def consistent(mode, tab):
    """every coloured vertex carries the id of an object of that colour's class"""
    cls_of = {(1 if r.list == "candidate" else -1) * int(r.obj_id): r.cls for _, r in tab.iterrows()}
    for h in ("lh", "rh"):
        lab, ctab, names = read_annot(f"{B}/review/PH6_surf_{mode}.{h}.annot")
        ids = nib.load(f"{B}/review/PH6_surf_ids.{h}.mgh").get_fdata().ravel()
        names = [n.decode() for n in names]
        for v in np.nonzero(lab >= 0)[0]:
            if cls_of.get(int(ids[v])) != names[lab[v]]: fail(f"{mode} {h}: vertex {v} colour/id mismatch")

# 1. colour by detector list
print("   " + run(f"{S}/surface_overlay.py", "PH6").stdout.splitlines()[1])
L = pd.read_csv(f"{B}/review/PH6_surf_list_objects.csv")
if (L.n_vertices == 0).any(): fail("objects not on the surface")
consistent("list", L)
want = {1: "S_k05", 2: "S_k07", 3: "S_k13", 6: "S_k09"}          # PH3 tram-track 1-3 and single bank 6
for t, sul in want.items():
    c = int(np.bincount(C[tr == t][C[tr == t] > 0]).argmax())
    got = L[(L.list == "candidate") & (L.obj_id == c)].sulcus.iloc[0]
    if got != sul: fail(f"cSS {t}: sulcus {got}, expected {sul}")
print(f"   list mode: annotation and id overlay consistent; cSS sulci {', '.join(want.values())} as constructed")

# 2. colour by expert calls (truth-derived answers: C cSS, VS surface vein / sinus, VC other)
run(f"{S}/expert_sheet.py", "PH6", "--include-excluded")
key = pd.read_csv(f"{B}/review/PH6_expert_key.csv"); ans = []
for _, r in key.iterrows():
    m = (C == r.obj_id) if r.list == "candidate" else (X == r.obj_id)
    ans.append(f"{r.sheet_no}{'C' if (tr[m] > 0).mean() > 0.2 else 'VS' if (sv[m] > 0).mean() > 0.5 else 'VC'}")
run(f"{S}/expert_import.py", "PH6", "expertS", "--letters", " ".join(ans))
print("   " + run(f"{S}/surface_overlay.py", "PH6", "--color-by", "expert").stdout.splitlines()[1])
E = pd.read_csv(f"{B}/review/PH6_surf_expert_objects.csv")
if (E.cls == "not_labelled").any(): fail("expert labels not matched to every object")
if set(E.cls) != {"cSS", "surface_vein", "other_vein"}: fail(f"expert classes {set(E.cls)}")
if E[E.cls == "cSS"].n_hidden.sum(): fail("a cSS call is hidden under another class")
consistent("expert", E)
r_ = run(f"{S}/surface_overlay.py", "PH6", "--color-by", "expert", "--reader", "nobody", ok=False)
if r_.returncode == 0 or "no labels from reader 'nobody'" not in r_.stderr: fail("unknown reader not refused clearly")
print(f"   expert mode: {E.cls.value_counts().to_dict()}; no cSS vertex hidden; colours consistent")
print("RESULT overlay OK")
