"""v4.17 (audit Stage 2) checks: recon-all label convention (sulcal CSF = 0), Hessian on the mask-filled image,
E8 artifact-zone rule only with T1 anatomy, stale FLAIR / Destrieux ignored, per-case manifest.
usage: check_stage2.py     (needs PH3 in $CSS_BASE)"""
import os, sys, json, subprocess, re, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
B = os.environ["CSS_BASE"]; S = f"{B}/scripts"; py = sys.executable
def fail(m): print("FAIL:", m); sys.exit(1)
def run(*a):
    r = subprocess.run([py, *a], capture_output=True, text=True)
    if r.returncode: print(r.stdout[-800:], r.stderr[-1500:]); fail("step failed: " + " ".join(a))
    return r.stdout
img = nib.load(f"{B}/data/PH3_swi.nii"); I = img.get_fdata().astype(np.float32); aff = img.affine
seg0 = nib.load(f"{B}/synthseg/PH3_seg.nii.gz").get_fdata().astype(np.int32)
brain = I > 0
def case(name, seg, t1=True, flair=None, data=None):
    nib.save(nib.Nifti1Image(I if data is None else data, aff), f"{B}/data/{name}_swi.nii")
    nib.save(nib.Nifti1Image(seg, aff), f"{B}/synthseg/{name}_seg.nii.gz")
    p = f"{B}/work/{name}_t1seg_swispace.nii.gz"
    if t1: nib.save(nib.Nifti1Image(seg, aff), p)
    elif os.path.exists(p): os.remove(p)
    if flair is not None: nib.save(nib.Nifti1Image(flair, aff), f"{B}/work/{name}_flair_swispace.nii.gz")
    run(f"{S}/align_seg.py", name); out = run(f"{S}/detect_css.py", name)
    return out, pd.read_csv(f"{B}/review/{name}_candidates.csv")
# FLAIR: CSF bright, tissue darker (no acute SAH anywhere)
tis = ndi.binary_fill_holes(brain & (seg0 > 0) & (seg0 != 24))
FL = np.where(brain & ~tis, 600.0, 300.0).astype(np.float32) * brain + \
     np.random.default_rng(0).normal(0, 10, I.shape).astype(np.float32) * brain

# 1. recon-all aparc+aseg convention: sulcal/extracerebral CSF unlabelled (0), label 24 only in a small basal region
seg_r = seg0.copy(); seg_r[seg_r == 24] = 0
basal = np.zeros_like(brain); basal[:, :, :12] = True
seg_r[(seg0 == 24) & basal] = 24
_, A = case("PH7A", seg0, flair=FL)               # SynthSeg convention
_, Bc = case("PH7B", seg_r, flair=FL)             # recon-all convention
key = ["vox_i", "vox_j", "vox_k", "volume_mm3"]
m = A.merge(Bc, on=key, suffixes=("_a", "_b"))
if len(m) != len(A) or len(m) != len(Bc): fail(f"candidates differ between label conventions ({len(A)} vs {len(Bc)})")
if (m.surface_contact_a - m.surface_contact_b).abs().max() > 0.01: fail("surface_contact depends on the CSF label convention")
if Bc.flair_csf_z.isna().all(): fail("FLAIR features switched off with recon-all labels")
if (m.flair_csf_z_a - m.flair_csf_z_b).abs().max() > 0.5: fail("FLAIR CSF reference depends on the label convention")
print(f"   label convention: {len(m)} identical candidates, surface_contact and FLAIR CSF z the same with sulcal CSF = 0")

# 2. Hessian on the mask-filled image (a sharp mask edge must not create curvature)
src = open(f"{S}/detect_css.py").read()
if src.count("gaussian_filter(If, sg, order=order)") != 2 or "gaussian_filter(Ic, sg, order=order)" in src:
    fail("3-D Hessian is not computed on the mask-filled image")
flat = np.where(brain, 200.0, 0.0); filled = np.where(brain, 200.0, 200.0)
sg = [1.6 / v for v in img.header.get_zooms()[:3]]
near_edge = brain & (ndi.distance_transform_edt(brain, sampling=img.header.get_zooms()[:3]) < 3)
e0 = np.abs(ndi.gaussian_filter(flat, sg, order=(2, 0, 0))[near_edge]).max()
e1 = np.abs(ndi.gaussian_filter(filled, sg, order=(2, 0, 0))[near_edge]).max()
if not e1 < 1e-2 * e0: fail("filled image still curves at the mask edge")   # residue = kernel truncation
print(f"   Hessian: second derivative at the mask edge {e0:.1f} (zeros outside) -> {e1:.3f} (filled)")

# 3. E8 artifact-zone rule: never on SWI-only anatomy
lut = f"{os.environ.get('FREESURFER_HOME', '')}/FreeSurferColorLUT.txt"
if os.path.exists(lut) and "lateralorbitofrontal" not in open(lut).read():
    open(lut, "a").write("1012 ctx-lh-lateralorbitofrontal 35 75 50 0\n2012 ctx-rh-lateralorbitofrontal 35 75 50 0\n")
seg_o = seg0.copy(); seg_o[seg0 == 1028] = 1012; seg_o[seg0 == 2028] = 2012
# a dark curvilinear line on the low frontal surface: relative height ~0.25 (above the 0.2 skull-base cut)
x, y, zz = np.meshgrid(*[(np.arange(n) - n / 2) * v for n, v in zip(I.shape, img.header.get_zooms()[:3])], indexing="ij")
dep = (1 - np.sqrt((x / 52) ** 2 + (y / 58) ** 2 + (zz / 44) ** 2)) * 52
ang = np.degrees(np.arctan2(y, x))
low = (dep > 0) & (dep < 1.2) & (zz > -24) & (zz < -18) & (np.abs(ang - 90) < 15) & (seg_o >= 1000)
I_o = I * np.where(low, 0.35, 1.0).astype(np.float32)
o1, _ = case("PH7C", seg_o, t1=False, data=I_o)
x1 = pd.read_csv(f"{B}/review/PH7C_excluded.csv")
if len(x1) and x1.reason.str.contains("artifact zone").any(): fail("E8 fired on SWI-only anatomy")
o2, _ = case("PH7D", seg_o, t1=True, data=I_o)
x2 = pd.read_csv(f"{B}/review/PH7D_excluded.csv")
n8 = int(x2.reason.str.contains("artifact zone").sum()) if len(x2) else 0
if n8 == 0: fail("E8 never fires even with T1 anatomy - the gating test is not exercising the rule")
print(f"   E8 artifact zone: 0 exclusions SWI-only, {n8} with T1 anatomy")

# 4. stale inputs: a FLAIR / Destrieux older than the SWI are ignored
import time
t0 = time.time()                                         # explicit times: file systems can share an mtime tick
os.utime(f"{B}/work/PH7A_flair_swispace.nii.gz", (t0 - 100, t0 - 100))   # FLAIR from before the SWI re-import
os.utime(f"{B}/data/PH7A_swi.nii", (t0 - 50, t0 - 50))
os.utime(f"{B}/work/PH7A_t1seg_swispace.nii.gz", (t0, t0))               # keep T1 labels current for this check
run(f"{S}/align_seg.py", "PH7A"); out = run(f"{S}/detect_css.py", "PH7A")
if "FLAIR IGNORED" not in out or not pd.read_csv(f"{B}/review/PH7A_candidates.csv").flair_csf_z.isna().all():
    fail("a stale FLAIR was used")
nib.save(nib.load(f"{B}/work/PH3_a2009s_swispace.nii.gz"), f"{B}/work/PH7A_a2009s_swispace.nii.gz")
os.utime(f"{B}/work/PH7A_a2009s_swispace.nii.gz", (t0 - 100, t0 - 100))
d = pd.read_csv(f"{B}/review/PH7A_candidates.csv"); d["accept"] = 1; d.to_csv(f"{B}/review/PH7A_candidates.csv", index=False)
sc = run(f"{S}/score_css.py", "PH7A")
if "Destrieux labels are older than the SWI" not in sc or "euclidean" not in sc: fail("stale Destrieux labels used")
print("   stale FLAIR and Destrieux labels ignored with a warning")

# 5. manifest
mf = json.load(open(f"{B}/review/PH7B_manifest.json"))["detect"]
need = {"pipeline_version", "inputs", "environment", "parameters", "qc_status"}
if not need <= set(mf): fail(f"manifest misses {need - set(mf)}")
import hashlib
if mf["inputs"]["swi"]["sha256"] != hashlib.sha256(open(f"{B}/data/PH7B_swi.nii", "rb").read()).hexdigest():
    fail("manifest SWI hash wrong")
if not {"swi", "labels", "t1_labels", "flair"} <= set(mf["inputs"]): fail(f"manifest inputs {list(mf['inputs'])}")
print(f"   manifest: {mf['pipeline_version']}, inputs {sorted(mf['inputs'])}, numpy {mf['environment']['numpy']}")
print("RESULT stage2 OK")
