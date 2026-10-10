"""v4.14 checks: case QC, generation-gate log, extent at the volume edge, CSV column contract.
usage: check_qc.py        (needs PH3 from phantom_v4.py in $CSS_BASE and the scripts in $CSS_BASE/scripts)
Each check prints one line and the script exits 1 on the first failure."""
import os, sys, json, subprocess, numpy as np, nibabel as nib, pandas as pd
B = os.environ["CSS_BASE"]; S = f"{B}/scripts"; py = sys.executable

V413 = ["cand_id", "hemi", "region", "label", "volume_mm3", "n_slices", "elongation", "surface_contact",
        "darkness_z", "artifact_zone", "midline_zone", "score", "vox_i", "vox_j", "vox_k", "cortex_frac",
        "cortex_dist_mm", "branch_per10mm", "surface_gradient", "vein_like", "long_structure", "pial_dist_mm",
        "bank_frac", "tube_ratio", "surface_alignment", "tram_frac", "vein_tree_mm", "mirror_dark_frac",
        "ich_dist_mm", "near_ich_suggest", "infratentorial", "flair_csf_z", "parenchyma_frac", "extent_mm",
        "cmb_like", "flair_bright", "flair_ctx_z", "flair_ctx_bright", "css_evidence", "vein_evidence",
        "edge_contrast", "infra_frac", "basal_frac", "rel_height", "vessel_ext_mm", "vessel_wm_mm",
        "vessel_run_mm", "axis_normal", "score_v3", "score_v4", "accept"]

def fail(m): print("FAIL:", m); sys.exit(1)
def run(*a):
    r = subprocess.run([py, *a], capture_output=True, text=True)
    if r.returncode: print(r.stdout[-800:], r.stderr[-1500:]); fail("step failed: " + " ".join(a))
    return r.stdout
def make_case(name, I, seg, aff, t1=True):
    """write a derived phantom; t1=True pretends the labels are T1-derived (fresh t1seg file)"""
    nib.save(nib.Nifti1Image(I.astype(np.float32), aff), f"{B}/data/{name}_swi.nii")
    nib.save(nib.Nifti1Image(seg.astype(np.int32), aff), f"{B}/synthseg/{name}_seg.nii.gz")
    p = f"{B}/work/{name}_t1seg_swispace.nii.gz"
    if t1: nib.save(nib.Nifti1Image(seg.astype(np.int32), aff), p)
    elif os.path.exists(p): os.remove(p)
    run(f"{S}/align_seg.py", name); run(f"{S}/detect_css.py", name)
    return json.load(open(f"{B}/review/{name}_qc.json"))

img = nib.load(f"{B}/data/PH3_swi.nii"); I0 = img.get_fdata().astype(np.float32); aff = img.affine
seg0 = nib.load(f"{B}/synthseg/PH3_seg.nii.gz").get_fdata().astype(np.int32)
truth = nib.load(f"{B}/work/PH3_truth.nii.gz").get_fdata().astype(int)

# 1. CSV column contract: v4.13 columns unchanged and first, new ones appended
run(f"{S}/detect_css.py", "PH3")
c = list(pd.read_csv(f"{B}/review/PH3_candidates.csv").columns)
if c[:len(V413)] != V413: fail(f"candidate columns changed: {c[:len(V413)]}")
if c[len(V413):] != ["qc_status", "rule_version", "n_pieces"]: fail(f"new columns: {c[len(V413):]}")
x = pd.read_csv(f"{B}/review/PH3_excluded.csv")
if list(x.columns[:2]) != ["excl_id", "reason"] or x.columns[-1] != "exclusion_confidence": fail("excluded columns")
if len(x) and not x.rule_version.str.match(r"E\d@v4\.\d").all(): fail("excluded rows without rule id")
print(f"columns: v4.13 order kept, appended {c[len(V413):]}; excluded rows carry rule ids {sorted(x.rule_version.unique())}")

# 2. generation-gate log: every dropped component has an allowed reason and never overlaps a shown/excluded one
d = pd.read_csv(f"{B}/work/PH3_dropped.csv")
dm = nib.load(f"{B}/work/PH3_dropped.nii.gz").get_fdata().astype(int)
cm = nib.load(f"{B}/review/PH3_candidates.nii.gz").get_fdata().astype(int)
xm = nib.load(f"{B}/review/PH3_excluded.nii.gz").get_fdata().astype(int)
ok_r = ("volume <", "<2 slices", "elongation <", "no cortex label")
if not len(d) or not d.reason.map(lambda r: r.startswith(ok_r)).all(): fail(f"dropped reasons {d.reason.unique()}")
if set(np.unique(dm)) - {0} != set(d.drop_id): fail("dropped map ids != dropped.csv ids")
if ((dm > 0) & ((cm > 0) | (xm > 0))).any(): fail("dropped voxels overlap candidates/excluded")
print(f"dropped log: {len(d)} components, reasons {d.reason.value_counts().to_dict()}")

# 3. orientation: same phantom stored with the S-I axis first -> insufficient_quality
perm = (2, 0, 1)
aff_t = aff.copy(); aff_t[:, :3] = aff[:, list(perm)]       # new axis i = old axis perm[i]
qt = make_case("PHT", np.transpose(I0, perm), np.transpose(seg0, perm), aff_t)
if qt["orientation_ok"] or qt["status"] != "insufficient_quality": fail(f"transposed phantom QC {qt['status']}")
print(f"orientation: transposed phantom -> {qt['status']} ({qt['reasons'][0]})")

# 4. label misregistration: labels shifted by 0 / 0.8 / 1.6 / 3.2 mm (along x, y and z) -> edge agreement grows
vals = []
for n in (0, 1, 2, 4):
    sh = (n, n, max(n * 2 // 3, 0))
    q = make_case(f"PHS{n}", I0, np.roll(seg0, sh, axis=(0, 1, 2)), aff)
    vals.append((n, q["edge_agreement_mm"], q["status"]))
e = [v[1] for v in vals]
print("label shift (voxels) -> edge agreement mm, status: " + ", ".join(f"{n}: {a} {s}" for n, a, s in vals))
if vals[0][2] != "ok": fail("unshifted T1-labelled phantom should be QC ok")
if not all(a < b for a, b in zip(e, e[1:])): fail(f"edge agreement not monotonic {e}")
if vals[-1][2] == "ok": fail("3.2 mm label shift still QC ok")

# 5. extent at the volume edge (B1): crop so that a truth lesion starts on array slice 0
k0 = int(np.argwhere(truth == 1)[:, 2].min()) + 1
aff_c = aff.copy(); aff_c[:3, 3] += aff[:3, 2] * k0
make_case("PHC", I0[:, :, k0:], seg0[:, :, k0:], aff_c)
cc = pd.read_csv(f"{B}/review/PHC_candidates.csv"); cmc = nib.load(f"{B}/review/PHC_candidates.nii.gz").get_fdata().astype(int)
xc = pd.read_csv(f"{B}/review/PHC_excluded.csv"); xmc = nib.load(f"{B}/review/PHC_excluded.nii.gz").get_fdata().astype(int)
vox = np.array(img.header.get_zooms()[:3]); edge_hits = 0
for tab, m, idc in ((cc, cmc, "cand_id"), (xc, xmc, "excl_id")):
    for _, r in tab.iterrows():
        p = np.argwhere(m == r[idc])
        if not len(p): continue
        ext = float(np.sqrt((((p.max(0) - p.min(0) + 1) * vox) ** 2).sum()))
        edge_hits += int(p[:, 2].min() == 0)
        if abs(ext - r.extent_mm) > 0.15: fail(f"extent_mm {r.extent_mm} != bounding box {ext:.1f} ({idc} {r[idc]})")
if not edge_hits: fail("no candidate touches slice 0 - test does not exercise the volume edge")
print(f"extent at the volume edge: {edge_hits} object(s) on slice 0, all extent_mm equal their bounding box")
print("RESULT qc OK")
