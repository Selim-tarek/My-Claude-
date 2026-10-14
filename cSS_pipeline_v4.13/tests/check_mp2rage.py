"""v4.18 mp2rage_robust.py: synthetic MP2RAGE (signed INV1, complex noise) -> background noise removed, brain kept.
usage: check_mp2rage.py    (scratch files in $CSS_BASE/work)"""
import os, sys, subprocess, numpy as np, nibabel as nib
B = os.environ["CSS_BASE"]; W = f"{B}/work"; os.makedirs(W, exist_ok=True)
def fail(m): print("FAIL:", m); sys.exit(1)
rng = np.random.default_rng(0); sh = (100, 110, 90)
x, y, z = np.meshgrid(*[np.arange(n) for n in sh], indexing="ij")
e = lambda a, b, c: ((x - 50) / a) ** 2 + ((y - 55) / b) ** 2 + ((z - 45) / c) ** 2 < 1
head, brain, wm = e(40, 45, 38), e(32, 37, 30), e(20, 24, 18)
i1 = np.where(head, -200.0, 0.0); i1[brain] = -80; i1[wm] = 150          # INV1 sign flips with T1
i2 = np.where(head, 600.0, 0.0); i2[brain] = 700; i2[wm] = 800
n = lambda v: v + rng.normal(0, 20, sh) + 1j * rng.normal(0, 20, sh)
c1, c2 = n(i1), n(i2)
uni = np.real(c1 * np.conj(c2)) / (np.abs(c1) ** 2 + np.abs(c2) ** 2)
pad = np.zeros(sh, bool); pad[8:-8, 8:-8, 8:-8] = True       # LEMON-like zero-padded / defaced border (v4.18.2)
for name, v in (("uni", np.round((uni + 0.5) * 4095)), ("inv1", np.abs(c1)), ("inv2", np.abs(c2))):
    nib.save(nib.Nifti1Image((v * pad).astype(np.float32), np.eye(4)), f"{W}/MP2_{name}.nii.gz")
r = subprocess.run([sys.executable, f"{B}/scripts/mp2rage_robust.py", *(f"{W}/MP2_{k}.nii.gz" for k in ("uni", "inv1", "inv2")),
                    f"{W}/MP2_robust.nii.gz"], capture_output=True, text=True)
if r.returncode: print(r.stdout, r.stderr); fail("mp2rage_robust failed")
U = nib.load(f"{W}/MP2_uni.nii.gz").get_fdata(); O = nib.load(f"{W}/MP2_robust.nii.gz").get_fdata()
from scipy import ndimage as ndi
band = (ndi.distance_transform_edt(~head) > 5) & (ndi.distance_transform_edt(~head) <= 15) & pad
hm = np.median(U[wm])
bg_before, bg_after = (U[band] > 0.25 * hm).mean(), (O[band] > 0.25 * hm).mean()
cor = float(np.corrcoef(U[brain], O[brain])[0, 1])
if not (bg_after < 0.05 and bg_before > 0.3 and cor > 0.95 and abs(np.median(O[wm]) - np.median(U[wm])) < 0.10 * np.median(U[wm])):
    fail(f"speckle {bg_before:.2f} -> {bg_after:.3f}, brain r {cor:.3f}, wm {np.median(O[wm]):.0f} vs {np.median(U[wm]):.0f}; {r.stdout.strip()}")
print(f"   robust MP2RAGE (zero-padded border): speckle 5-15 mm outside the head {bg_before:.0%} -> {bg_after:.1%}, "
      f"brain correlation {cor:.3f}, white matter within 10 %; {r.stdout.strip().split('; ')[0].split(': ')[-1]}")
print("RESULT mp2rage OK")
