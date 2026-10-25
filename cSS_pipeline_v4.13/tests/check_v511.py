"""v5.1.1 checks: per-piece evidence of merged objects, CSS_MIXED_REVIEW, valid-sample fractions, missing-input
indicators, versioned expert labels.
usage: check_v511.py    (unit parts need nothing; the rest needs PH3 from phantom_v4.py, aligned, in $CSS_BASE)
1. css_pieces: a cSS-like band merged with a vein-like tube -> mixed_object 1; one kind only, or a vein piece under
   15 % of the volume -> 0; no merging -> no piece evidence. rescue_mixed only for E5 / E6 / E7, only when enabled.
2. PH3: the new columns are present; CSS_MIXED_REVIEW=1 changes nothing when no object is mixed.
3. a lopsided head mask (part of one hemisphere cut away): mirror points outside the brain are reported by
   mirror_valid_frac < 1 and are not counted as dark; ref_valid_frac within 0..1.
4. missing-input columns match the inputs (PH3: SWI-only anatomy, no FLAIR, no phase); the ranker's add_missing.
5. expert labels: label_schema_version, label_class, confidence_missing written; --require-confidence refuses.
PHANTOM / UNIT CHECKS ONLY."""
import os, sys, json, shutil, subprocess, numpy as np, nibabel as nib, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import css_pieces
from train_ranker import add_missing
def fail(m): print("FAIL:", m); sys.exit(1)

# 1. units
sh = (20, 20, 10); mm = np.zeros(sh, bool); pieces = np.zeros(sh, int)
pial = np.zeros(sh, np.float32); tube = np.full(sh, np.nan, np.float32)
mm[2:14, 2:7, 4] = True; pieces[2:14, 2:7, 4] = 1; pial[2:14, 2:7, 4] = -0.5; tube[2:14, 2:7, 4] = 0.2   # band 60 vox
mm[14:18, 9:11, 3:8] = True; pieces[14:18, 9:11, 3:8] = 2; pial[14:18, 9:11, 3:8] = 1.0; tube[14:18, 9:11, 3:8] = 0.8  # tube 40
e = css_pieces.piece_evidence(mm, pieces, pial, tube)
if e["mixed_object"] != 1 or abs(e["piece_vein_frac"] - 0.4) > 0.01: fail(f"band + tube: {e}")
only = mm & (pieces == 1)
if css_pieces.piece_evidence(only, pieces, pial, tube)["mixed_object"] != 0: fail("one kind of piece called mixed")
small = mm.copy(); small[14:18, 9:11, 4:8] = False                   # tube piece 8 of 68 voxels = 12 %
if css_pieces.piece_evidence(small, pieces, pial, tube)["mixed_object"] != 0: fail("12 % vein piece called mixed")
nm = css_pieces.piece_evidence(mm, None, pial, tube)
if nm["mixed_object"] != 0 or not np.isnan(nm["piece_vein_frac"]): fail(f"no merging: {nm}")
R = css_pieces.rescue_mixed
if not (R("E5@v4.8", 1, True) and R("E7@v4.12", 1, True) and not R("E1@v4.9.2", 1, True)
        and not R("E5@v4.8", 0, True) and not R("E5@v4.8", 1, False)):
    fail("rescue_mixed logic")
print("pieces: band + tube -> mixed (vein fraction 0.40); one kind / 12 % vein / no merging -> not mixed; "
      "rescue only for E5-E7 with CSS_MIXED_REVIEW=1")
d = add_missing(pd.DataFrame(dict(phase_z=[np.nan, 1.0], flair_csf_z=[np.nan, np.nan], flair_ctx_z=[np.nan, 2.0],
                                  surf_width_mm=[np.nan, 3.0])))
if d.missing_phase_features.tolist() != [1, 0] or d.missing_flair_features.tolist() != [1, 0] \
        or d.missing_surface_features.tolist() != [1, 0]:
    fail(f"add_missing {d.to_dict('list')}")
print("ranker missing-input indicators: phase / FLAIR / surface families flagged per row before imputation")

B = os.environ.get("CSS_BASE")
if not B or not os.path.exists(f"{B}/data/PH3_swi.nii"):
    print("PH3 parts skipped (no PH3 in $CSS_BASE)"); print("RESULT v511 OK"); sys.exit(0)
S = f"{B}/scripts"; py = sys.executable
def run(*a, env=None, ok=True):
    e_ = dict(os.environ); e_.pop("CSS_MIXED_REVIEW", None); e_.update(env or {})
    r = subprocess.run([py, *a], capture_output=True, text=True, env=e_)
    if ok and r.returncode: print(r.stdout[-800:], r.stderr[-1200:]); fail("step failed: " + " ".join(a))
    return r
# 2. columns, MIXED_REVIEW no-op
run(f"{S}/detect_css.py", "PH3"); c0, x0 = (pd.read_csv(f"{B}/review/PH3_{f}.csv") for f in ("candidates", "excluded"))
new = ["piece_vol_vox", "piece_pial_mm", "piece_tube", "piece_vein_frac", "mixed_object", "mixed_object_review",
       "merge_gap_mm", "ref_valid_frac", "mirror_valid_frac", "region_n_vox", "missing_t1_anatomy",
       "missing_flair_features", "missing_phase_features"]
if any(c not in c0.columns or c not in x0.columns for c in new): fail("v5.1.1 columns missing")
run(f"{S}/detect_css.py", "PH3", env={"CSS_MIXED_REVIEW": "1"}); c1 = pd.read_csv(f"{B}/review/PH3_candidates.csv")
if int(c0.mixed_object.sum()) == 0 and not c0[["vox_i", "vox_j", "vox_k", "volume_mm3"]].equals(c1[["vox_i", "vox_j", "vox_k", "volume_mm3"]]):
    fail("CSS_MIXED_REVIEW changed the output although no object is mixed")
print(f"PH3: v5.1.1 columns present in candidates and excluded; {int(c0.mixed_object.sum())} mixed objects; "
      f"CSS_MIXED_REVIEW=1 changes nothing here")
# 4. missing inputs on PH3 (SWI-only, no FLAIR, no phase)
if not (c0.missing_t1_anatomy.eq(1).all() and c0.missing_flair_features.eq(1).all() and c0.missing_phase_features.eq(1).all()):
    fail("PH3 missing-input columns wrong")
# 3. lopsided head mask
img = nib.load(f"{B}/data/PH3_swi.nii"); I = img.get_fdata().astype(np.float32); aff = img.affine
seg = nib.load(f"{B}/synthseg/PH3_seg.nii.gz").get_fdata().astype(np.int32)
I2, s2 = I.copy(), seg.copy(); cut = int(0.72 * I.shape[0]); I2[cut:] = 0; s2[cut:] = 0     # part of one side gone
nib.save(nib.Nifti1Image(I2, aff), f"{B}/data/PHA_swi.nii"); nib.save(nib.Nifti1Image(s2, aff), f"{B}/synthseg/PHA_seg.nii.gz")
p = f"{B}/work/PHA_t1seg_swispace.nii.gz"
if os.path.exists(p): os.remove(p)
run(f"{S}/align_seg.py", "PHA"); run(f"{S}/detect_css.py", "PHA")
a = pd.concat([pd.read_csv(f"{B}/review/PHA_{f}.csv") for f in ("candidates", "excluded")])
low = a[a.mirror_valid_frac < 1]
if not len(low): fail("no candidate with mirror points outside the cut brain")
if ((low.mirror_valid_frac < 0.5) & low.mirror_dark_frac.notna()).any(): fail("mirror dark computed from mostly outside points")
if not a.ref_valid_frac.between(0, 1).all(): fail("ref_valid_frac outside 0..1")
print(f"lopsided head mask: {len(low)} objects with mirror points outside the brain (mirror_valid_frac "
      f"{low.mirror_valid_frac.min():.2f}-{low.mirror_valid_frac.max():.2f}); below 0.5 -> mirror_dark_frac NaN, not 'dark'")
# 5. expert labels
for f in ("candidates", "excluded"):
    shutil.copy(f"{B}/review/PH3_{f}.csv", f"{B}/review/PHK_{f}.csv")
c = pd.read_csv(f"{B}/review/PHK_candidates.csv").head(3)
key = c[["vox_i", "vox_j", "vox_k", "volume_mm3"]].assign(sheet_no=[1, 2, 3], list="candidate", obj_id=c.cand_id.values, cand_id=c.cand_id.values)
key.to_csv(f"{B}/review/PHK_expert_key.csv", index=False)
r = run(f"{S}/expert_import.py", "PHK", "rdrA", "--letters", "1C:4 2VS 3N:5", "--require-confidence", ok=False)
if r.returncode == 0: fail("--require-confidence accepted an answer without confidence")
run(f"{S}/expert_import.py", "PHK", "rdrA", "--letters", "1C:4 2VS 3N:5")
L = pd.read_csv(f"{B}/labels/candidate_labels.csv"); L = L[(L.study_id == "PHK") & (L.reader_id == "rdrA")]
if not {"label_schema_version", "label_class", "confidence_missing"} <= set(L.columns): fail("label schema columns missing")
if L.label_schema_version.iloc[0] != "css-labels/1.0" or L.sort_values("sheet_no").confidence_missing.tolist() != [0, 1, 0]:
    fail(f"label rows {L[['sheet_no', 'label_schema_version', 'confidence_missing']].to_dict('list')}")
if L[L.sheet_no == 2].label_class.iloc[0] != "surface_vein": fail("label_class not written")
print("expert labels: css-labels/1.0, fine label + class + confidence_missing written; --require-confidence refuses")
print("RESULT v511 OK")
