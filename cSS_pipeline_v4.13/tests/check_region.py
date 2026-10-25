"""v4.19.2 checks: regional + contralateral reference (detect_css.py region_z, region_dark_frac, mirror_region_z,
cand_mirror_z). usage: check_region.py   (needs PH3 from phantom_v4.py, aligned, in $CSS_BASE)
PH3 + one gyrus darkened by 25 % and given its own label (1024 / 2024, like a naturally dark motor strip in
recon-all's aparc): PHR = both hemispheres (physiological, symmetric), PHU = left only.
PHANTOM REGRESSION CHECKS ONLY: they do not validate any threshold for real data."""
import os, sys, subprocess, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
B = os.environ["CSS_BASE"]; S = f"{B}/scripts"; py = sys.executable
def fail(m): print("FAIL:", m); sys.exit(1)
img = nib.load(f"{B}/data/PH3_swi.nii"); I0 = img.get_fdata().astype(np.float32); aff = img.affine
seg0 = nib.load(f"{B}/synthseg/PH3_seg.nii.gz").get_fdata().astype(np.int32)
a2 = nib.load(f"{B}/work/PH3_a2009s_swispace.nii.gz").get_fdata().astype(int)
tr = nib.load(f"{B}/work/PH3_truth.nii.gz").get_fdata().astype(int)
GL = (a2 == 11104) & (seg0 >= 1000); GR = (a2 == 12104) & (seg0 >= 1000)   # gyrus sector 3, its own mirror; no lesions
if (tr[GL | GR] > 0).any() or GL.sum() < 1000: fail("test gyrus overlaps a lesion / is missing")
def case(name, dark):
    I = I0.copy(); seg = seg0.copy(); seg[GL] = 1024; seg[GR] = 2024
    I *= 1 - 0.25 * ndi.gaussian_filter(dark.astype(np.float32), 1.0)
    nib.save(nib.Nifti1Image(I, aff), f"{B}/data/{name}_swi.nii")
    for p in (f"{B}/synthseg/{name}_seg.nii.gz", f"{B}/work/{name}_t1seg_swispace.nii.gz"):
        nib.save(nib.Nifti1Image(seg, aff), p)
    for st in ("align_seg.py", "detect_css.py"):
        r = subprocess.run([py, f"{S}/{st}", name], capture_output=True, text=True)
        if r.returncode: print(r.stdout[-800:], r.stderr[-1500:]); fail(f"{st} {name} failed")
    rows = []
    for lst, idc in (("candidates", "cand_id"), ("excluded", "excl_id")):
        d = pd.read_csv(f"{B}/review/{name}_{lst}.csv"); M = nib.load(f"{B}/review/{name}_{lst}.nii.gz").get_fdata().astype(int)
        for _, r_ in d.iterrows():
            m = M == r_[idc]
            k = "cSS" if (tr[m] > 0).mean() > 0.2 else ("dark gyrus" if dark[m].mean() > 0.3 else "other")
            rows.append(dict(kind=k, **{c: r_[c] for c in ("region_z", "region_dark_frac", "mirror_region_z", "cand_mirror_z")}))
    D = pd.DataFrame(rows); c, g = D[D.kind == "cSS"], D[D.kind == "dark gyrus"]
    if len(c) < 4 or not len(g): fail(f"{name}: expected cSS and dark-gyrus objects, got {D.kind.value_counts().to_dict()}")
    print(f"{name}: cSS region_z {c.region_z.max():.2f}..{c.region_z.min():.2f}, dark frac <= {c.region_dark_frac.max():.3f}, "
          f"mirror {c.mirror_region_z.abs().max():.2f} | dark gyrus region_z {g.region_z.min():.2f}..{g.region_z.max():.2f}, "
          f"dark frac {g.region_dark_frac.min():.2f}..{g.region_dark_frac.max():.2f}, mirror {g.mirror_region_z.min():.2f}..{g.mirror_region_z.max():.2f}")
    if (c.region_z > -4).any() or (c.region_dark_frac > 0.1).any(): fail(f"{name}: cSS should be focal (far darker than its region)")
    if (g.region_z.abs() > 2).any() or (g.region_dark_frac < 0.3).any(): fail(f"{name}: dark gyrus should read as broadly dark")
    if (c.mirror_region_z.abs() > 1).any(): fail(f"{name}: unchanged regions should be symmetric")
    return g
g = case("PHR", GL | GR)
if (g.mirror_region_z.abs() > 1).any(): fail("bilateral dark gyrus should be symmetric")
g = case("PHU", GL)
if (g.mirror_region_z > -2).any(): fail("one-sided dark gyrus should be darker than its mirror region")
print("RESULT region OK")
