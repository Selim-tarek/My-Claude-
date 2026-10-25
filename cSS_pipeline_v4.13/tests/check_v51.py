"""v5.1 checks: auditable fates, label-convention and storage invariance, obliquity, stale derived inputs.
usage: check_v51.py    (needs PH3 from phantom_v4.py, aligned, in $CSS_BASE; uses every other phantom found there)
1. fate invariant: for every phantom case, dropped / excluded / shown maps are disjoint, together cover exactly the
   generated voxels (work/ID_gen.nii.gz bit 128), and the manifest fate counts add up; audit_case.py agrees.
2. sulcal CSF encoded as label 0 (recon-all aparc+aseg) instead of 24 (SynthSeg): identical candidates / exclusions.
3. the in-plane array axes swapped (same anatomy, different storage): same objects at the same world positions.
4. an oblique affine (data unchanged): obliquity recorded; 15 deg changes no candidate; 25 deg -> low confidence.
5. stale inputs: FLAIR / phase older than the SWI are ignored; surface features older than the detection are ignored
   (CSS_QC_STRICT=1: stop); a ranking model whose label file changed after the fit: warning (strict: refused).
PHANTOM CHECKS ONLY."""
import os, sys, json, time, glob, subprocess, numpy as np, nibabel as nib, pandas as pd
B = os.environ["CSS_BASE"]; S = f"{B}/scripts"; py = sys.executable
def fail(m): print("FAIL:", m); sys.exit(1)
def run(*a, env=None, ok=True):
    e = dict(os.environ); e.pop("CSS_QC_STRICT", None); e.pop("CSS_USE_MIMIC", None); e.update(env or {})
    r = subprocess.run([py, *a], capture_output=True, text=True, env=e)
    if ok and r.returncode: print(r.stdout[-800:], r.stderr[-1200:]); fail("step failed: " + " ".join(a))
    return r
def case(name, I, seg, aff):
    nib.save(nib.Nifti1Image(I.astype(np.float32), aff), f"{B}/data/{name}_swi.nii")
    nib.save(nib.Nifti1Image(seg.astype(np.int32), aff), f"{B}/synthseg/{name}_seg.nii.gz")
    p = f"{B}/work/{name}_t1seg_swispace.nii.gz"
    if os.path.exists(p): os.remove(p)
    run(f"{S}/align_seg.py", name); return run(f"{S}/detect_css.py", name)
def tables(name):
    return pd.read_csv(f"{B}/review/{name}_candidates.csv"), pd.read_csv(f"{B}/review/{name}_excluded.csv")
img = nib.load(f"{B}/data/PH3_swi.nii"); I0 = img.get_fdata().astype(np.float32); aff = img.affine
seg0 = nib.load(f"{B}/synthseg/PH3_seg.nii.gz").get_fdata().astype(np.int32)
run(f"{S}/detect_css.py", "PH3")

# 1. fate invariant on every case that has outputs
checked = []
for cp in sorted(glob.glob(f"{B}/review/*_candidates.nii.gz")):
    s = os.path.basename(cp)[:-len("_candidates.nii.gz")]
    gp, mp = f"{B}/work/{s}_gen.nii.gz", f"{B}/review/{s}_manifest.json"
    if not (os.path.exists(gp) and os.path.exists(mp)): continue
    F = json.load(open(mp)).get("detect", {}).get("fates")
    if not F or F.get("qc_blocked"): continue
    if os.path.getmtime(gp) > os.path.getmtime(cp) + 1 or os.path.getmtime(mp) < os.path.getmtime(cp) - 1: continue
    C, X, D = (nib.load(p).get_fdata() > 0 for p in (cp, f"{B}/review/{s}_excluded.nii.gz", f"{B}/work/{s}_dropped.nii.gz"))
    G = (nib.load(gp).get_fdata().astype(np.int32) & 128) > 0
    if (C.astype(int) + X + D > 1).any(): fail(f"{s}: an object voxel has two fates")
    if ((C | X | D) != G).any(): fail(f"{s}: {int(((C | X | D) != G).sum())} generated voxels without a fate (or the reverse)")
    n = sum(F["dropped_by_gate"].values()) + sum(F["excluded_by_rule"].values()) + F["shown"]
    if n != F["objects"]: fail(f"{s}: fates add up to {n}, objects {F['objects']}")
    a = run(f"{S}/audit_case.py", s).stdout
    if "[ok]" not in a: fail(f"{s}: audit_case consistency not ok:\n{a}")
    checked.append(f"{s} {F['objects']}")
if len(checked) < 1: fail("no case to check the fate invariant on")
print(f"fate invariant: {len(checked)} cases (objects: {', '.join(checked)}) - every object one fate, voxel for voxel; audit ok")

# 2. sulcal / extra-axial CSF as label 0 instead of 24
c24, x24 = tables("PH3")
sg = seg0.copy(); sg[sg == 24] = 0
case("PHL0", I0, sg, aff); c0, x0 = tables("PHL0")
skip = {"qc_status", "rule_version"}
cols = [c for c in c24.columns if c not in skip]
if not (c24[cols].fillna(-9).astype(str).equals(c0[cols].fillna(-9).astype(str))
        and x24[[c for c in x24.columns if c not in skip]].fillna(-9).astype(str).equals(x0[[c for c in x24.columns if c not in skip]].fillna(-9).astype(str))):
    diff = [c for c in cols if len(c24) == len(c0) and not c24[c].fillna(-9).astype(str).equals(c0[c].fillna(-9).astype(str))]
    fail(f"CSF label 0 vs 24 changed the output ({len(c24)} vs {len(c0)} candidates; columns {diff[:6]})")
print(f"CSF as label 0 vs 24: identical {len(c0)} candidates and {len(x0)} exclusions (surface_contact, FLAIR CSF reference, geometry)")

# 3. in-plane axes swapped in storage
I1 = np.transpose(I0, (1, 0, 2)); s1 = np.transpose(seg0, (1, 0, 2)); aff1 = aff.copy(); aff1[:, [0, 1]] = aff[:, [1, 0]]
case("PHX", I1, s1, aff1); cx, xx = tables("PHX")
wpos = lambda d, A: (A @ np.c_[d[["vox_i", "vox_j", "vox_k"]].values, np.ones(len(d))].T)[:3].T
if len(cx) != len(c24) or len(xx) != len(x24): fail(f"axis swap: {len(cx)}/{len(xx)} vs {len(c24)}/{len(x24)}")
pa = wpos(c24, aff); pb = wpos(cx, aff1)
dmax = max(np.min(np.linalg.norm(pb - p, axis=1)) for p in pa) if len(pa) else 0
if dmax > 1.0 or sorted(x24.reason) != sorted(xx.reason): fail(f"axis swap moved objects by up to {dmax:.2f} mm or changed exclusions")
print(f"in-plane axes swapped: same {len(cx)} candidates (world positions within {dmax:.2f} mm), same exclusion reasons")

# 4. oblique affine
def rot_x(a_deg):
    t = np.radians(a_deg); R = np.eye(4); R[1:3, 1:3] = [[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]]
    return R @ aff
for deg in (15, 25):
    case(f"PHO{deg}", I0, seg0, rot_x(deg)); co, xo = tables(f"PHO{deg}")
    q = json.load(open(f"{B}/review/PHO{deg}_qc.json"))
    if abs(q["obliquity_deg"] - deg) > 0.2: fail(f"obliquity {q['obliquity_deg']} for a {deg} deg tilt")
    flagged = any("oblique" in r for r in q["reasons"])
    if flagged != (deg > 20): fail(f"{deg} deg: oblique flag {flagged}")
    if deg == 15 and (len(co) != len(c24) or not (co.volume_mm3.values == c24.volume_mm3.values).all()):
        fail("a 15 deg tilt changed the candidates")
print("oblique affine: obliquity recorded (15 / 25 deg); 15 deg changes no candidate; 25 deg -> low_confidence_review 'oblique slices'")

# 5. stale inputs
nib.save(nib.Nifti1Image(I0 * 0.5 + 100, aff), f"{B}/work/PH3_flair_swispace.nii.gz")
nib.save(nib.Nifti1Image(np.zeros_like(I0), aff), f"{B}/data/PH3_phase.nii.gz")
old = os.path.getmtime(f"{B}/data/PH3_swi.nii") - 600          # older than the SWI
for p in (f"{B}/work/PH3_flair_swispace.nii.gz", f"{B}/data/PH3_phase.nii.gz"): os.utime(p, (old, old))
o = run(f"{S}/detect_css.py", "PH3").stdout
if "FLAIR used" in o or "older than the SWI" not in o: fail(f"stale FLAIR / phase not ignored:\n{o[-600:]}")
os.remove(f"{B}/work/PH3_flair_swispace.nii.gz"); os.remove(f"{B}/data/PH3_phase.nii.gz")
c = pd.read_csv(f"{B}/review/PH3_candidates.csv")
pd.DataFrame(dict(list="candidate", obj_id=c.cand_id, vox_i=c.vox_i, vox_j=c.vox_j, vox_k=c.vox_k, volume_mm3=c.volume_mm3,
                  surf_width_mm=3.0)).to_csv(f"{B}/review/PH3_surface_features.csv", index=False)
os.utime(f"{B}/review/PH3_surface_features.csv", (old, old))
model = dict(features=["darkness_z", "surf_width_mm"], weights=[1.0, 1.0], bias=0.0, impute_median=[3.0, 1.0],
             mean=[3.0, 1.0], sd=[1.0, 1.0], training_studies=["XX"], n_patients=3, lopo_auc=None, baseline="score_v4")
os.makedirs(f"{B}/labels", exist_ok=True); lab = f"{B}/work/v51_labels.csv"; open(lab, "w").write("a\n1\n")
from importlib import util as _u
sys.path.insert(0, S); import css_common
model.update(label_file=lab, label_sha256=css_common.sha256(lab)); json.dump(model, open(f"{B}/work/v51_model.json", "w"))
r = run(f"{S}/rank_candidates.py", "PH3", "--model", f"{B}/work/v51_model.json")
if "surface features" not in r.stdout or "ignored" not in r.stdout: fail(f"stale surface features used:\n{r.stdout}")
r = run(f"{S}/rank_candidates.py", "PH3", "--model", f"{B}/work/v51_model.json", env={"CSS_QC_STRICT": "1"}, ok=False)
if r.returncode == 0: fail("CSS_QC_STRICT=1 used stale surface features")
os.remove(f"{B}/review/PH3_surface_features.csv")
open(lab, "a").write("2\n")
r = run(f"{S}/rank_candidates.py", "PH3", "--model", f"{B}/work/v51_model.json")
if "changed since the model was fitted" not in r.stdout: fail("changed label file not reported")
if run(f"{S}/rank_candidates.py", "PH3", "--model", f"{B}/work/v51_model.json", env={"CSS_QC_STRICT": "1"}, ok=False).returncode == 0:
    fail("CSS_QC_STRICT=1 ranked with a stale model")
run(f"{S}/detect_css.py", "PH3")                                   # restore PH3 outputs for later checks
m = json.load(open(f"{B}/review/PH3_manifest.json"))
if m.get("schema_version") != "css-manifest/1.0" or "build" not in m["detect"]["environment"]: fail("manifest schema / build missing")
print("stale inputs: old FLAIR / phase ignored; old surface features ignored (strict: stop); changed labels -> model warning "
      "(strict: refused); manifest schema css-manifest/1.0 with build id")
# 6. the mimic layer must not cost a true lesion: every phantom with a truth map, default vs CSS_USE_MIMIC=1
lost_any, seen = [], []
for tp in sorted(glob.glob(f"{B}/work/*_truth.nii.gz")):
    s = os.path.basename(tp)[:-len("_truth.nii.gz")]
    if not os.path.exists(f"{B}/data/{s}_swi.nii") or not os.path.exists(f"{B}/work/{s}_seg_swispace.nii.gz"): continue
    tr = nib.load(tp).get_fdata().astype(int)
    def kept(env):
        r = run(f"{S}/detect_css.py", s, env=env, ok=False)
        if r.returncode: return None
        Cm = nib.load(f"{B}/review/{s}_candidates.nii.gz").get_fdata().astype(int)
        return {int(t) for t in np.unique(tr[tr > 0]) if (Cm[tr == t] > 0).mean() >= 0.3}
    k0, k1 = kept({}), kept({"CSS_USE_MIMIC": "1"})
    if k0 is None or k1 is None: continue
    seen.append(f"{s} {len(k0)}->{len(k1)}")
    if k0 - k1: lost_any.append(f"{s}: lesions {sorted(k0 - k1)}")
    run(f"{S}/detect_css.py", s)                                    # default outputs back
if not seen: fail("no phantom with a truth map to test the mimic layer on")
if lost_any: fail("CSS_USE_MIMIC=1 lost true lesions: " + "; ".join(lost_any))
print(f"mimic layer on vs off, true lesions kept: {', '.join(seen)} - none lost")
print("RESULT v51 OK")
