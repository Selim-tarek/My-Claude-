"""v5.2 checks: surface-aligned patches (surface_aligned_patch.py), the patch classifier (classify_surface_patches.py)
and the shared model code (css_model.py). Synthetic volumes and tables with known answers; no FreeSurfer, no phantoms.
usage: check_v52.py    (writes into $CSS_BASE/v52_check/)
1. frame: a dark SHEET under a flat pial surface vs a dark TUBE in the CSF above it. n points outward (tissue on -n,
   CSF on +n) for both; the sheet's dark voxels span the patch across the surface (v), the tube's do not; the sheet
   lies below the surface, the tube above. Right-handed (u, v, n).
2. a reversed signed-distance map is flipped by the labels (normal_flipped) and gives the same n; an oblique,
   anisotropic affine (30 deg about x, 1 x 1 x 2 mm) gives n = the rotated slice axis; < 3 voxels -> skipped.
3. CLI on a synthetic study: phase absent, FLAIR off the SWI grid, then FLAIR older than the SWI -> each recorded as
   missing with its reason, never zero-filled (valid all False, feature columns NaN, <ch>_missing 1); patches near
   the array edge are partial (0 < valid fraction < 1) with the fixed shape; voxels outside the brain mask are never
   valid; the 2-voxel object is listed as skipped with its reason (no silent loss).
4. css_model: a fold whose training part has one class is skipped and logged; out-of-fold scores of a patient do not
   change when that patient's own labels change (patient grouping, plain and nested).
5. classifier on synthetic feature tables (6 dev + 1 test study, phase missing in half): train fits on dev only,
   ablation A-G reported, model hash stored; --split test refused; too few studies refused; predict refuses a
   training study, scores the test study (cSS above veins), refuses a tampered model and a schema mismatch.
6. PH3 (when present in $CSS_BASE): patches for every shown and excluded object of the real detector output, one
   manifest row per object (none lost), features extracted."""
import os, sys, json, shutil, subprocess, numpy as np, nibabel as nib, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); S = os.path.join(HERE, "..", "scripts"); sys.path.insert(0, S)
import surface_aligned_patch as SAP, classify_surface_patches as CSP, css_model as CM
from css_common import cand_uid
py = sys.executable
B0 = os.environ.get("CSS_BASE", "/tmp"); B = f"{B0}/v52_check"
def fail(m): print("FAIL:", m); sys.exit(1)

# ---------- synthetic volume: flat pial surface at k = 31.5, tissue below, CSF above ----------
N = 64; ii, jj, kk = np.meshgrid(*[np.arange(N)] * 3, indexing="ij")
labels = np.where(kk < 32, 3, np.where(kk < 56, 24, 0)).astype(np.int16)
swi = np.where(kk < 56, 100.0, 0.0).astype(np.float32); swi[ii < 14] = 0.0        # i < 14: outside the brain
sheet = (ii >= 16) & (ii <= 30) & (jj >= 20) & (jj <= 32) & (kk >= 30) & (kk <= 31)
tube = (ii >= 16) & (ii <= 36) & (jj >= 45) & (jj <= 47) & (kk >= 33) & (kk <= 34)
tiny = (ii == 40) & (jj == 10) & (kk == 30) | (ii == 41) & (jj == 10) & (kk == 30)
edge = (ii >= 40) & (ii <= 52) & (jj <= 4) & (kk >= 30) & (kk <= 31)
for m in (sheet, tube, edge): swi[m] = 20.0
swi[tiny] = 20.0
psd = (kk - 31.5).astype(np.float32); tissue = (labels > 0) & (labels != 24)
A1 = np.eye(4)

def dark_extent(z, fr, axis):
    """extent (mm, 0.5 mm grid) of dark samples along one patch axis"""
    G = SAP.grid(fr, [24, 12, 8], [0.5] * 3); v, ins = SAP.sample(z, A1, G)
    d = ins & (v < 50); idx = np.nonzero(d.any(axis=tuple(a for a in range(3) if a != axis)))[0]
    return (idx.max() - idx.min()) * 0.5 if len(idx) else 0.0, G, d

fs = SAP.frame(sheet, psd, A1, np.ones(3), tissue); ft = SAP.frame(tube, psd, A1, np.ones(3), tissue)
for name, fr in (("sheet", fs), ("tube", ft)):
    if fr["n"] @ np.array([0, 0, 1.0]) < 0.99: fail(f"{name}: normal {fr['n']} is not outward (+z)")
    if abs(np.linalg.det(np.c_[fr["u"], fr["v"], fr["n"]]) - 1) > 1e-6: fail(f"{name}: frame not right-handed")
ev_s = dark_extent(swi, fs, 1)[0]; ev_t = dark_extent(swi, ft, 1)[0]
if not (ev_s >= 10 and ev_t <= 3.5): fail(f"across-surface dark extent sheet {ev_s} mm, tube {ev_t} mm")
if not (fs["pial_sd_at_center_mm"] < 0 < ft["pial_sd_at_center_mm"]): fail("sheet should sit below, tube above the surface")
lab_s = CSP.channel_features("labels", SAP.sample(labels.astype(np.float32), A1, SAP.grid(fs, [24, 12, 8], [0.5] * 3), 0)[0],
                             np.ones((49, 25, 17), bool))[1]
if not (lab_s[2] == 1.0 and lab_s[9] == 0.0): fail(f"labels along n: tissue fraction first/last bin {lab_s[2]}/{lab_s[9]}")
print(f"1. frame: n outward for sheet and tube, right-handed; dark extent across the surface sheet {ev_s:.1f} mm vs tube "
      f"{ev_t:.1f} mm; pial distance at centre sheet {fs['pial_sd_at_center_mm']:+.1f} / tube {ft['pial_sd_at_center_mm']:+.1f} mm; "
      f"tissue fraction along n {lab_s[2]:.0f} -> {lab_s[9]:.0f}")

# ---------- 2. sign inversion, oblique affine, < 3 voxels ----------
ff = SAP.frame(sheet, -psd, A1, np.ones(3), tissue)
if not ff["normal_flipped"] or ff["n"] @ fs["n"] < 0.99: fail(f"reversed distance map: flipped {ff['normal_flipped']}, n {ff['n']}")
c, s_ = np.cos(np.pi / 6), np.sin(np.pi / 6); R = np.array([[1, 0, 0], [0, c, -s_], [0, s_, c]])
Ao = np.eye(4); Ao[:3, :3] = R @ np.diag([1, 1, 2.0]); Ao[:3, 3] = [-30, -20, -60]
fo = SAP.frame(sheet, (kk - 31.5).astype(np.float32) * 2, Ao, np.array([1, 1, 2.0]), tissue)
if fo["n"] @ R[:, 2] < 0.99: fail(f"oblique affine: n {np.round(fo['n'], 3)} vs expected {np.round(R[:, 2], 3)}")
try:
    SAP.frame(tiny, psd, A1, np.ones(3), tissue); fail("2-voxel object was not refused")
except ValueError as e:
    msg = str(e)
print(f"2. reversed map flipped by the labels (n . n_ref {ff['n'] @ fs['n']:.3f}); oblique 1x1x2 mm affine: n . rotated "
      f"slice axis {fo['n'] @ R[:, 2]:.3f}; tiny object refused ('{msg}')")

# ---------- 3. CLI on a synthetic study ----------
shutil.rmtree(B, ignore_errors=True)
for d in ("data", "work", "review", "labels"): os.makedirs(f"{B}/{d}")
def save(a, p, aff=A1): nib.save(nib.Nifti1Image(a, aff), p)
save(swi, f"{B}/data/T1S_swi.nii")
save(psd, f"{B}/work/T1S_pial_sd.nii.gz"); save(labels, f"{B}/work/T1S_seg_swispace.nii.gz")
Aoff = A1.copy(); Aoff[0, 3] = 3.0
save(swi, f"{B}/work/T1S_flair_swispace.nii.gz", Aoff)                   # FLAIR not on the SWI grid
objs = [sheet, tube, tiny, edge]; M = np.zeros(swi.shape, np.int16)
for n_, m in enumerate(objs, 1): M[m] = n_
save(M, f"{B}/review/T1S_candidates.nii.gz")
rows = [dict(cand_id=n_, vox_i=int(np.argwhere(m).mean(0)[0]), vox_j=int(np.argwhere(m).mean(0)[1]),
             vox_k=int(np.argwhere(m).mean(0)[2]), volume_mm3=float(m.sum())) for n_, m in enumerate(objs, 1)]
pd.DataFrame(rows).to_csv(f"{B}/review/T1S_candidates.csv", index=False)
env = dict(os.environ, CSS_BASE=B)
def run(*a, ok=True):
    r = subprocess.run([py, *a], capture_output=True, text=True, env=env)
    if ok and r.returncode: print(r.stdout[-1500:], r.stderr[-1500:]); fail("step failed: " + " ".join(a[:2]))
    return r
run(f"{S}/surface_aligned_patch.py", "T1S", "--lists", "candidates")
man = pd.read_csv(f"{B}/review/T1S_patches.csv")
if len(man) != 4: fail(f"manifest has {len(man)} rows for 4 objects (silent loss)")
sk = man[man.status == "skipped"]
if list(sk.obj_id) != [3] or "< 3" not in str(sk.skip_reason.iloc[0]): fail(f"skipped rows {sk.to_dict('records')}")
if "flair:not on the SWI grid" not in str(man.missing_channels.iloc[0]) or "phase:absent" not in str(man.missing_channels.iloc[0]):
    fail(f"missing channels {man.missing_channels.iloc[0]}")
z = np.load(f"{B}/review/T1S_patches/candidates_1.npz")
if "image_phase" in z.files or "image_flair" in z.files or z["valid_phase"].any() or z["valid_flair"].any():
    fail("a missing channel was written as an image / valid voxels")
W = SAP.to_voxel(A1, z["world_grid_mm"])
if (z["valid_swi"] & (W[..., 0] < 13.5)).any(): fail("voxels outside the brain mask are valid")
if not (z["valid_swi"].any() and (~z["valid_swi"] & (W[..., 0] < 13.5)).any()): fail("sheet patch should reach outside the brain")
ze = np.load(f"{B}/review/T1S_patches/candidates_4.npz")
if not (0 < ze["valid_swi"].mean() < 1) or ze["valid_swi"].shape != z["valid_swi"].shape: fail("edge patch should be partial")
run(f"{S}/classify_surface_patches.py", "features", "T1S")
F = pd.read_csv(f"{B}/review/T1S_patch_features.csv")
if not ((F.phase_missing == 1).all() and (F.flair_missing == 1).all() and (F.swi_missing == 0).all()
        and F.phase_mean.isna().all() and F.flair_p50.isna().all()): fail("missing channels not encoded as missing + NaN")
# FLAIR on the grid but older than the SWI -> stale, also missing
save(swi, f"{B}/work/T1S_flair_swispace.nii.gz"); t = os.path.getmtime(f"{B}/data/T1S_swi.nii") - 600
os.utime(f"{B}/work/T1S_flair_swispace.nii.gz", (t, t))
run(f"{S}/surface_aligned_patch.py", "T1S", "--lists", "candidates")
man2 = pd.read_csv(f"{B}/review/T1S_patches.csv")
if "flair:older than the SWI" not in str(man2.missing_channels.iloc[0]): fail(f"stale FLAIR: {man2.missing_channels.iloc[0]}")
print(f"3. CLI: 4 objects -> {int((man.status == 'ok').sum())} patches + 1 skipped ('{sk.skip_reason.iloc[0]}'); phase absent, "
      f"FLAIR off-grid then stale -> missing with reason, NaN features, <ch>_missing 1; outside-brain voxels never valid; "
      f"edge patch valid fraction {ze['valid_swi'].mean():.2f} with the fixed shape {ze['valid_swi'].shape}")

# ---------- 4. css_model folds ----------
rng = np.random.default_rng(3)
g = np.repeat(["a", "b", "c"], 10); y = np.r_[np.tile([0, 1], 5), np.zeros(20, int)]; X = rng.normal(size=(30, 2)) + y[:, None]
sc, log = CM.lopo(X, y, g)
st = {f["study_id"]: f["status"] for f in log}
if st["a"] != "skipped" or not np.isnan(sc[g == "a"]).all(): fail(f"one-class training fold not skipped: {log}")
g2 = np.repeat(list("pqrst"), 12); y2 = np.tile([0, 0, 1], 20); X2 = rng.normal(size=(60, 3)) + y2[:, None]
for nested in (False, True):
    s1 = CM.lopo(X2, y2, g2, nested=nested)[0]; y3 = y2.copy(); y3[g2 == "p"] = 1 - y3[g2 == "p"]
    s2 = CM.lopo(X2, y3, g2, nested=nested)[0]
    if not np.allclose(s1[g2 == "p"], s2[g2 == "p"]): fail(f"patient p's own labels changed its out-of-fold scores (nested {nested})")
print("4. css_model: one-class training fold skipped and logged; a patient's out-of-fold scores ignore its own labels "
      "(plain and nested)")

# ---------- 5. classifier on synthetic feature tables ----------
names = []
for ch in CSP.CHANNELS: names += CSP.channel_features(ch, None, None)[0]
meta, lab = [], []
for p in range(7):
    s = f"Q{p:02d}"; n = 20; yy = (np.arange(n) % 3 == 0).astype(int)
    d = pd.DataFrame(rng.normal(size=(n, len(names))), columns=names)
    for ch in CSP.CHANNELS: d[f"{ch}_missing"] = 0.0; d[f"{ch}_valid_fraction"] = 0.9
    d["swi_mean"] = -2.5 * yy + rng.normal(0, 0.3, n)                  # cSS darker on average (strong: a plumbing test)
    if p % 2:                                                           # phase missing in half the studies
        d[[c for c in names if c.startswith("phase_")]] = np.nan; d["phase_missing"] = 1.0; d["phase_valid_fraction"] = 0.0
    r = pd.DataFrame(dict(vox_i=np.arange(n), vox_j=p, vox_k=1, volume_mm3=30.0))
    d.insert(0, "obj_id", np.arange(1, n + 1)); d.insert(0, "list", "candidates")
    d.insert(0, "cand_uid", [cand_uid(s, rr) for _, rr in r.iterrows()]); d.insert(0, "study_id", s)
    d["feature_schema"] = CSP.FEATURE_SCHEMA; d.to_csv(f"{B}/review/{s}_patch_features.csv", index=False)
    meta.append(dict(study_id=s, split="test" if p == 6 else "dev"))
    lab += [dict(study_id=s, cand_uid=u, consensus_class="cSS" if t else "surface_vein", method="agree") for u, t in zip(d.cand_uid, yy)]
pd.DataFrame(meta).to_csv(f"{B}/labels/study_meta.csv", index=False); pd.DataFrame(lab).to_csv(f"{B}/labels/consensus.csv", index=False)
if run(f"{S}/classify_surface_patches.py", "train", "--split", "test", ok=False).returncode == 0: fail("--split test accepted")
if run(f"{S}/classify_surface_patches.py", "train", "--min-studies", "10", ok=False).returncode == 0: fail("too few studies accepted")
out = run(f"{S}/classify_surface_patches.py", "train", "--no-nested").stdout
Mo = json.load(open(f"{B}/review/surface_patch_model.json"))
if "Q06" in Mo["training_studies"] or len(Mo["training_studies"]) != 6: fail(f"training studies {Mo['training_studies']}")
if set(Mo["ablations"]) != set(CSP.ABLATIONS): fail("ablations missing")
aG = Mo["ablations"]["G_all"]["auc"]; aA = Mo["ablations"]["A_swi"]["auc"]
if not (aA > 0.85 and aG > 0.75): fail(f"LOPO AUC A_swi {aA}, G_all {aG}")
if Mo["model_hash"] != CM.model_hash(Mo) or not Mo["research_only"]: fail("model hash / research_only")
r = run(f"{S}/classify_surface_patches.py", "predict", "Q00", "Q06").stdout
if "used to train the model - not predicted" not in r or os.path.exists(f"{B}/review/Q00_patch_predictions.csv"): fail("training study predicted")
P = pd.read_csv(f"{B}/review/Q06_patch_predictions.csv").merge(pd.DataFrame(lab), on=["study_id", "cand_uid"])
a6 = CM.auc((P.consensus_class == "cSS").astype(int), P.css_score)
if a6 < 0.8 or not (P.research_only == 1).all() or P.model_hash.iloc[0] != Mo["model_hash"]: fail(f"test-study AUC {a6} / provenance")
Mt = dict(Mo); Mt["weights"] = [w + 0.1 for w in Mo["weights"]]; json.dump(Mt, open(f"{B}/review/tampered.json", "w"))
if run(f"{S}/classify_surface_patches.py", "predict", "Q06", "--model", f"{B}/review/tampered.json", ok=False).returncode == 0:
    fail("tampered model accepted")
Ms = dict(Mo); Ms["feature_schema"] = "css-patch-features/0.9"; Ms["model_hash"] = CM.model_hash(Ms)
json.dump(Ms, open(f"{B}/review/oldschema.json", "w"))
if run(f"{S}/classify_surface_patches.py", "predict", "Q06", "--model", f"{B}/review/oldschema.json", ok=False).returncode == 0:
    fail("schema mismatch accepted")
print(f"5. classifier: fitted on 6 dev studies (test study left out), {len(Mo['feature_names'])} features; LOPO AUC A_swi {aA:.3f}, "
      f"G_all {aG:.3f}; held-out study AUC {a6:.3f}; --split test, too few studies, a training study, a tampered model and a "
      f"schema mismatch all refused")

PB = os.environ.get("CSS_BASE")
if PB and os.path.exists(f"{PB}/work/PH3_seg_swispace.nii.gz"):
    e2 = dict(os.environ)
    r = subprocess.run([py, f"{S}/detect_css.py", "PH3"], capture_output=True, text=True, env=e2)
    if r.returncode: fail("detect_css PH3")
    r = subprocess.run([py, f"{S}/surface_aligned_patch.py", "PH3"], capture_output=True, text=True, env=e2)
    if r.returncode: print(r.stdout[-800:], r.stderr[-800:]); fail("surface_aligned_patch PH3")
    n_obj = sum(len(pd.read_csv(f"{PB}/review/PH3_{l}.csv")) for l in ("candidates", "excluded"))
    man = pd.read_csv(f"{PB}/review/PH3_patches.csv")
    if len(man) != n_obj: fail(f"PH3: {n_obj} objects but {len(man)} manifest rows")
    r = subprocess.run([py, f"{S}/classify_surface_patches.py", "features", "PH3"], capture_output=True, text=True, env=e2)
    if r.returncode: fail("patch features PH3")
    print(f"6. PH3: {n_obj} objects -> {int((man.status == 'ok').sum())} patches, {int((man.status != 'ok').sum())} skipped with a "
          f"reason; {r.stdout.strip()}")
else:
    print("6. phantom part skipped")
print("RESULT v52 OK")
