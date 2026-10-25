"""v4.18 physics-based synthetic lesions: dipole field (css_common.susceptibility_signal) against analytic values,
then make_synthetic.py --physics on PH1: blooming grows with TE, calcium has the opposite phase, lesions detected.
usage: check_physics.py    (needs PH1 from phantom.py, aligned, in $CSS_BASE)"""
import os, sys, re, shutil, subprocess, numpy as np, nibabel as nib, pandas as pd
B = os.environ["CSS_BASE"]; S = f"{B}/scripts"; py = sys.executable
sys.path.insert(0, S)
from css_common import susceptibility_signal, GAMMA_HZ_T
def fail(m): print("FAIL:", m); sys.exit(1)
def run(*a):
    r = subprocess.run([py, *a], capture_output=True, text=True)
    if r.returncode: print(r.stdout[-600:], r.stderr[-1500:]); fail("step failed: " + " ".join(a))
    return r.stdout

# 1. analytic checks: long cylinder inside field = chi (3 cos^2 th - 1) / 6; sphere outside: + along B0, - across
sh, vox, aff = (80, 80, 80), (0.5, 0.5, 0.5), np.diag([0.5, 0.5, 0.5, 1.0])
k = 2 * np.pi * GAMMA_HZ_T * 3.0 * 0.02 * 1e-6
for name, sl, th in (("parallel", (slice(39, 41), slice(39, 41), slice(5, 75)), 0.0),
                     ("perpendicular", (slice(5, 75), slice(39, 41), slice(39, 41)), np.pi / 2)):
    chi = np.zeros(sh); chi[sl] = 0.45
    ph, _ = susceptibility_signal(chi, vox, aff, 3.0, 0.02)
    want = k * 0.45 * (3 * np.cos(th) ** 2 - 1) / 6; got = float(np.median(ph[sl]))
    if abs(got - want) > 0.1 * abs(want): fail(f"cylinder {name}: phase {got:.2f} vs analytic {want:.2f}")
x, y, z = np.meshgrid(*[np.arange(80)] * 3, indexing="ij")
chi = np.zeros(sh); chi[(x - 40) ** 2 + (y - 40) ** 2 + (z - 40) ** 2 < 9] = 1.0
ph, _ = susceptibility_signal(chi, vox, aff, 3.0, 0.02)
if not (ph[40, 40, 46] > 0 > ph[46, 40, 40] and abs(ph[40, 40, 46] / ph[46, 40, 40] + 2) < 0.4):
    fail(f"sphere dipole pattern {ph[40, 40, 46]:.2f} / {ph[46, 40, 40]:.2f} (expect ratio -2)")
# B0 follows the affine: the same cylinder stored with axes swapped (array axis 0 = scanner z) reads parallel
aff2 = np.zeros((4, 4)); aff2[0, 2] = aff2[1, 1] = aff2[2, 0] = 0.5; aff2[3, 3] = 1
chi = np.zeros(sh); chi[5:75, 39:41, 39:41] = 0.45
ph, _ = susceptibility_signal(chi, vox, aff2, 3.0, 0.02)
if abs(float(np.median(ph[5:75, 39:41, 39:41])) - k * 0.45 / 3) > 0.1 * k * 0.45 / 3: fail("B0 direction not taken from the affine")
print(f"   dipole model: cylinders match the analytic field (within 10 %), sphere pattern +2:-1, B0 from the affine")

# 2. make_synthetic --physics on PH1: TE dependence, calcium phase sign, detection
shutil.copy(f"{B}/synthseg/PH1_seg.nii.gz", f"{B}/synthseg/PH1S_seg.nii.gz")
I = nib.load(f"{B}/data/PH1_swi.nii").get_fdata(); res = {}
for te in (10, 40):
    run(f"{S}/make_synthetic.py", "PH1", "8", "0.6", "1", "--physics", "--te-ms", str(te), "--calcium", "2", "--veins", "4")
    Sx = nib.load(f"{B}/data/PH1S_swi.nii").get_fdata(); t = nib.load(f"{B}/work/PH1S_truth.nii.gz").get_fdata() > 0
    v = nib.load(f"{B}/work/PH1S_veins.nii.gz").get_fdata().astype(int); ph = nib.load(f"{B}/work/PH1S_dphase.nii.gz").get_fdata()
    kd = pd.read_csv(f"{B}/work/PH1S_vein_kinds.csv").set_index("vein_id").kind
    ca = np.isin(v, [i for i, x in kd.items() if x == "calcium"])
    if not ca.any(): fail("no calcium decoys inserted")
    res[te] = dict(bloom=int(((Sx < 0.8 * I) & (I > 0) & ~t & (v == 0)).sum()), ca_phase=float(np.median(ph[ca])),
                   lesion=float(np.median(Sx[t] / I[t])))
if not res[40]["bloom"] > 1.5 * res[10]["bloom"]: fail(f"blooming does not grow with TE: {res}")
if not (res[10]["ca_phase"] < 0 and res[40]["lesion"] < res[10]["lesion"]): fail(f"calcium phase / TE darkening wrong: {res}")
run(f"{S}/make_synthetic.py", "PH1", "8", "0.6", "1", "--physics", "--calcium", "2", "--veins", "4", "--surface-veins", "2")
run(f"{S}/align_seg.py", "PH1S"); run(f"{S}/detect_css.py", "PH1S")
r = re.search(r"found=(\d+) total=(\d+)", run(f"{S}/eval_synthetic.py", "PH1S"))
if int(r.group(1)) < 6: fail(f"physics lesions: only {r.group(1)}/{r.group(2)} found")
print(f"   physics lesions on PH1: darkened voxels outside objects TE 10 ms {res[10]['bloom']} -> TE 40 ms {res[40]['bloom']}; "
      f"calcium phase {res[10]['ca_phase']:+.2f} rad (diamagnetic); found {r.group(1)}/{r.group(2)} at TE 20 ms")
print("RESULT physics OK")
