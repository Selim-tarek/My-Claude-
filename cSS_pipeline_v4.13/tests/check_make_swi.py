"""v4.18 make_swi.py: synthetic GRE (wrapped background phase, paramagnetic veins, diamagnetic calcium) with
both phase sign conventions -> sign found automatically, veins darkened, calcium untouched, grid kept.
usage: check_make_swi.py    (uses $CSS_BASE/work as scratch)"""
import os, sys, subprocess, numpy as np, nibabel as nib
from scipy import ndimage as ndi
B = os.environ["CSS_BASE"]; W = f"{B}/work"; os.makedirs(W, exist_ok=True)
def fail(m): print("FAIL:", m); sys.exit(1)
sh = (120, 140, 90); vox = 0.8
x, y, z = np.meshgrid(*[np.arange(n) for n in sh], indexing="ij")
head = ((x - 60) / 50) ** 2 + ((y - 70) / 60) ** 2 + ((z - 45) / 38) ** 2 < 1
rng = np.random.default_rng(0); veins = np.zeros(sh, bool)
for k in range(15):
    p0 = np.array([rng.integers(25, 95), rng.integers(25, 115), rng.integers(20, 70)]); d = rng.normal(size=3); d /= np.linalg.norm(d)
    pts = (p0 + np.outer(np.linspace(0, 25, 200), d)).round().astype(int)
    ok = np.all((pts >= 0) & (pts < sh), 1); veins[tuple(pts[ok].T)] = True
veins = ndi.binary_dilation(veins) & head
calc = (x - 60) ** 2 + (y - 40) ** 2 + (z - 45) ** 2 < 9
M = np.where(head, 1000.0, 0.0); M[veins] *= 0.6; M[calc] *= 0.6
M = np.clip(M + rng.normal(0, 20, sh) * head, 0, None)
aff = np.diag([-vox, vox, -vox, 1.0]); aff[:3, 3] = [50, -50, 40]          # LAI, like LEMON
for s in (1, -1):
    bg = 0.08 * (x - 60) + 0.05 * (y - 70) + 0.002 * (z - 45) ** 2          # wraps several times
    ph = np.angle(np.exp(1j * (bg + s * 1.2 * ndi.gaussian_filter(veins.astype(float), 0.7) / 0.35 - s * 1.2 * calc))) * head
    nib.save(nib.Nifti1Image(M.astype(np.float32), aff), f"{W}/GRE_mag.nii.gz")
    nib.save(nib.Nifti1Image(((ph + np.pi) / (2 * np.pi) * 4095 - 2048).astype(np.float32), aff), f"{W}/GRE_pha.nii.gz")
    r = subprocess.run([sys.executable, f"{B}/scripts/make_swi.py", f"{W}/GRE_mag.nii.gz", f"{W}/GRE_pha.nii.gz",
                        f"{W}/GRE_swi.nii.gz"], capture_output=True, text=True)
    if r.returncode: print(r.stdout, r.stderr); fail("make_swi failed")
    if f"paramagnetic sign {s:+d}" not in r.stdout: fail(f"phase sign {s:+d} not found: {r.stdout}")
    S = nib.load(f"{W}/GRE_swi.nii.gz")
    if S.shape != sh or not np.allclose(S.affine, aff): fail("grid / affine changed")
    S = S.get_fdata(); rv = float(np.median(S[veins] / M[veins].clip(1))); rc = float(np.median(S[calc] / M[calc].clip(1)))
    if not (rv < 0.3 and rc > 0.95): fail(f"sign {s:+d}: veins x{rv:.2f}, calcium x{rc:.2f}")
    print(f"   sign {s:+d}: found; veins x{rv:.2f} of magnitude, calcium x{rc:.2f}")
print("RESULT make_swi OK")
