"""v5.0 checks: mimic layer (scripts/mimic_rules.py) and its hooks in detect_css.py.
usage: check_mimic.py    (unit part needs nothing; the integration part needs PH3 from phantom_v4.py, aligned, in $CSS_BASE)
1. phase rule: a PARAMAGNETIC sheet / tube / blob in any orientation to B0 gives phase_z > 0, a DIAMAGNETIC one < 0,
   for weak and phase-wrapping susceptibilities (dipole forward model css_common.susceptibility_signal + noise).
   Also shows why the sign of the median phase is NOT used: it is negative for a paramagnetic sheet flat to B0.
2. scanner sign convention fitted on paramagnetic veins, for both conventions.
3. microbleed-cluster rule: 4 compact parenchymal blobs -> cluster; thin cSS-like bands and a single blob -> not.
4. FLAIR rule classes.
5. PHM = PH3 + a physics phase image + one cSS lesion made DIAMAGNETIC (calcified band) + a row of 4 microbleeds:
   CSS_USE_MIMIC unset -> same candidates and exclusions as without the phase (recorded only);
   CSS_USE_MIMIC=1 -> the calcified band is excluded as E10, every paramagnetic cSS lesion is kept.
PHANTOM CHECKS ONLY: they show the rules measure what they are meant to on constructed data, not that the
thresholds hold on real scans."""
import os, sys, subprocess, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import mimic_rules as M
from css_common import susceptibility_signal, homodyne_highpass
def fail(m): print("FAIL:", m); sys.exit(1)
rng = np.random.default_rng(2)
iso = np.array([0.8, 0.8, 0.8]); aff_iso = np.diag([0.8, 0.8, 0.8, 1.0])

def rejector(hp, shape, conv=1, vox=iso, aff=aff_iso):
    r = M.MimicRejector.__new__(M.MimicRejector)
    r.vox, r.si, r.aff, r.brain, r.hp, r.conv = np.asarray(vox, float), 2, aff, np.ones(shape, bool), hp, conv
    return r

# ---- 1. phase rule across orientations, strengths, shapes
n = 96; zz, yy, xx = np.mgrid[:n, :n, :n]
def obj(kind):
    m = np.zeros((n, n, n), bool)
    if kind == "sheet flat": m[30:54, 30:54, 48] = True
    elif kind == "sheet vertical": m[48, 30:54, 30:54] = True
    elif kind == "tube || B0": m[47:49, 47:49, 30:66] = True
    elif kind == "tube _|_ B0": m[30:66, 47:49, 47:49] = True
    else: m[((xx - 48) ** 2 + (yy - 48) ** 2 + (zz - 48) ** 2) <= 4] = True
    return m
full = tuple(slice(0, n) for _ in range(3))
lines, naive_flat = [], None
for chi in (0.15, -0.15, 0.5, -0.5):
    for kind in ("sheet flat", "sheet vertical", "tube || B0", "tube _|_ B0", "blob"):
        m = obj(kind)
        ph, _ = susceptibility_signal(m * chi, iso, aff_iso, 3.0, 0.02, pad=24)
        ph = np.angle(np.exp(1j * (ph + rng.normal(0, 0.15, ph.shape))))
        hp = homodyne_highpass(ph, iso, 2, 4.0)
        pz, pr = rejector(hp, m.shape).phase_evidence(m, full)
        if not np.isfinite(pz) or (pz > 0) != (chi > 0):
            fail(f"phase rule: chi {chi:+} {kind} gave phase_z {pz} (wrong sign)")
        if chi == 0.5 and kind == "sheet flat": naive_flat = float(np.median(hp[m]))
        lines.append(f"{chi:+.2f} {kind}: {pz:+.1f}")
print("phase rule, 20 objects (chi ppm, shape: phase_z): " + "; ".join(lines))
ph, _ = susceptibility_signal(obj("sheet flat") * 0.5, iso, aff_iso, 3.0, 0.02, pad=24)
med = float(np.median(homodyne_highpass(ph, iso, 2, 4.0)[obj("sheet flat")]))
ph0, _ = susceptibility_signal(obj("sheet flat") * 0.15, iso, aff_iso, 3.0, 0.02, pad=24)
med0 = float(np.median(ph0[obj("sheet flat")]))
if not med0 < 0: fail("expected the raw phase inside a paramagnetic sheet flat to B0 to be negative")
print(f"why not the median sign: paramagnetic sheet flat to B0 (+0.15 ppm) has median phase {med0:+.2f} rad (negative); "
      f"at +0.5 ppm the filtered median is {med:+.2f} rad (wrapped)")

# ---- 2. sign convention fitted on veins
ves = np.zeros((n, n, n), bool)
for _ in range(8):
    c = rng.integers(20, n - 20, 3); d = rng.normal(size=3); d /= np.linalg.norm(d)
    for t in np.linspace(-14, 14, 70):
        p = np.round(c + t * d).astype(int); ves[p[0] - 1:p[0] + 1, p[1] - 1:p[1] + 1, p[2] - 1:p[2] + 1] = True
ph, _ = susceptibility_signal(ves * 0.45, iso, aff_iso, 3.0, 0.02, pad=24)
for sign in (1, -1):
    hp = homodyne_highpass(np.angle(np.exp(1j * (sign * ph + rng.normal(0, 0.15, ph.shape)))), iso, 2, 4.0)
    pred = M.highpass_linear(M.dipole_field(ves.astype(np.float32), iso, aff_iso), iso, 2)
    a, t, r = M.fit_amplitude(hp[ves], pred[ves])
    if np.sign(a) != sign or abs(t) < M.CONV_MIN_T: fail(f"sign convention {sign:+d} fitted as {np.sign(a):+.0f} (t {t:.1f})")
    print(f"scanner convention {sign:+d}: fitted {int(np.sign(a)):+d}, t {t:.1f}, r {r:.2f}")

# ---- 3. microbleed-cluster rule, at 0.8 iso, 0.8 x 0.8 x 1.2 and 0.5 x 0.5 x 2 mm
for vx in ([0.8, 0.8, 0.8], [0.8, 0.8, 1.2], [0.5, 0.5, 2.0]):
    vx = np.array(vx); sh = (90, 90, 40); Pm = np.indices(sh).transpose(1, 2, 3, 0) * vx; T = np.ones(sh, bool)
    c = np.array(sh) * vx / 2; row = np.zeros(sh, bool)
    for k in range(4): row |= ((Pm - (c + [0, (k - 1.5) * 5, 0])) ** 2).sum(-1) <= 1.6 ** 2
    s_, n_, _ = M.mb_cluster_rule(row, T, vx)
    if n_ != 4 or s_ < 0.5: fail(f"{vx}: 4 microbleeds (3.2 mm) counted as {n_} (score {s_})")
    bands = []
    for th in (1, 2, 3, 4):
        flat = np.zeros(sh, bool); flat[10:80, 10:80, 18:18 + th] = True
        vert = np.zeros(sh, bool); vert[10:80, 44:44 + th, 5:35] = True
        bands += [M.mb_cluster_rule(flat, T, vx)[1], M.mb_cluster_rule(vert, T, vx)[1]]
    one = ((Pm - c) ** 2).sum(-1) <= 1.6 ** 2
    if max(bands) >= 3 or M.mb_cluster_rule(one, T, vx)[1] >= 3: fail(f"{vx}: band or single blob called a cluster {bands}")
    print(f"microbleed cluster {vx.tolist()} mm: 4 blobs -> {n_} (score {s_}); cSS-like bands 1-4 voxels thick -> max {max(bands)}")

# ---- 4. FLAIR rule
cases = {"none": (np.nan, np.nan), "acute": (4.0, 1.0), "thrombosis": (1.0, 4.5), "mixed": (5.0, 5.0)}
for want, (fz, fcz) in cases.items():
    sc, cls, _ = M.flair_rule(fz, fcz)
    if cls != want or (want != "none") != (sc > 0): fail(f"FLAIR rule {fz}, {fcz} -> {cls} {sc}, expected {want}")
print("FLAIR rule: none / acute (CSF z 4) / thrombosis (cortex z 4.5) / mixed - as expected; z <= 3 gives none")

# ---- 5. integration on PHM (PH3 + phase + calcified band + microbleeds)
B = os.environ.get("CSS_BASE")
if not B or not os.path.exists(f"{B}/data/PH3_swi.nii"):
    print("integration part skipped (no PH3 in $CSS_BASE)"); print("RESULT mimic OK"); sys.exit(0)
S = f"{B}/scripts"; py = sys.executable
img = nib.load(f"{B}/data/PH3_swi.nii"); I = img.get_fdata().astype(np.float32); aff = img.affine
vox = np.array(img.header.get_zooms()[:3]); seg = nib.load(f"{B}/synthseg/PH3_seg.nii.gz").get_fdata().astype(np.int32)
tr = nib.load(f"{B}/work/PH3_truth.nii.gz").get_fdata().astype(int)
vn = nib.load(f"{B}/work/PH3_veins.nii.gz").get_fdata().astype(int)
CALC = 5                                                         # convexity crown lesion -> calcified band (4 is never generated on PH3)
chi = np.zeros(I.shape, np.float32); chi[(tr > 0) & (tr != CALC)] = 0.3; chi[tr == CALC] = -0.3; chi[vn > 0] = 0.45
# 4 microbleeds (r 1.6 mm, 5 mm apart: separate blobs that the 3 mm merge gap may join into one candidate)
# at the cortex / white-matter junction (lobar microbleeds are cortico-subcortical; >= half surrounded by parenchyma)
ctx = (seg >= 1000) & (tr == 0) & (vn == 0) & ~ndi.binary_dilation(tr > 0, iterations=6)
csf_d = ndi.distance_transform_edt(~np.isin(seg, [0, 24]), sampling=vox)
cand = np.argwhere(ctx & ndi.binary_dilation(np.isin(seg, [2, 41]), iterations=1) & (csf_d >= 2.5))
c0 = cand[len(cand) // 3].astype(float)
gz = np.gradient(ndi.gaussian_filter((seg >= 1000).astype(np.float32), 2))
tdir = np.cross([g[tuple(c0.astype(int))] for g in gz], [0, 0, 1.0]); tdir /= np.linalg.norm(tdir) + 1e-9
P = np.argwhere(np.ones(I.shape, bool)).reshape(I.shape + (3,)) * vox
mb = np.zeros(I.shape, bool)
for k in range(4):
    cc = (c0 + tdir * k * 5.0 / vox) * vox
    mb |= ((P - cc) ** 2).sum(-1) <= 1.6 ** 2
chi[mb] = 0.8
I2 = I.copy(); I2[mb] *= 0.35
ph, att = susceptibility_signal(chi, vox, aff, 3.0, 0.02, pad=16)
I2 = I2 * att                                                    # blooming, as on a real GRE
ph = np.angle(np.exp(1j * (ph + rng.normal(0, 0.1, ph.shape)))).astype(np.float32)
for name in ("PHM0", "PHM"):
    nib.save(nib.Nifti1Image(I2, aff), f"{B}/data/{name}_swi.nii")
    nib.save(nib.Nifti1Image(seg, aff), f"{B}/synthseg/{name}_seg.nii.gz")
    subprocess.run([py, f"{S}/align_seg.py", name], check=True, capture_output=True)
nib.save(nib.Nifti1Image(ph, aff), f"{B}/data/PHM_phase.nii.gz")      # PHM0: same scan without a phase image
if os.path.exists(f"{B}/data/PHM0_phase.nii.gz"): os.remove(f"{B}/data/PHM0_phase.nii.gz")
def det(name, use):
    env = dict(os.environ); env.pop("CSS_USE_MIMIC", None)
    if use: env["CSS_USE_MIMIC"] = "1"
    r = subprocess.run([py, f"{S}/detect_css.py", name], capture_output=True, text=True, env=env)
    if r.returncode: print(r.stdout[-800:], r.stderr[-1500:]); fail(f"detect_css {name} failed")
    out = {}
    for lst, idc in (("candidates", "cand_id"), ("excluded", "excl_id")):
        d = pd.read_csv(f"{B}/review/{name}_{lst}.csv"); Mk = nib.load(f"{B}/review/{name}_{lst}.nii.gz").get_fdata().astype(int)
        d["lesion"] = [int(np.bincount(tr[Mk == v][tr[Mk == v] > 0]).argmax()) if (tr[Mk == v] > 0).mean() > 0.2 else 0
                       for v in d[idc]]
        d["mb"] = [bool(mb[Mk == v].mean() > 0.3) for v in d[idc]]
        out[lst] = d
    return out, r.stdout
a0, _ = det("PHM0", False); a1, o1 = det("PHM", False); a2, o2 = det("PHM", True)
key = ["vox_i", "vox_j", "vox_k", "volume_mm3"]
if not (a0["candidates"][key].equals(a1["candidates"][key]) and a0["excluded"][key + ["reason"]].equals(a1["excluded"][key + ["reason"]])):
    fail("with CSS_USE_MIMIC unset, adding a phase image changed candidates or exclusions")
c1 = a1["candidates"]; calc = pd.concat([c1[c1.lesion == CALC], a1["excluded"][a1["excluded"].lesion == CALC]])
para = c1[(c1.lesion > 0) & (c1.lesion != CALC)]
print("phase_z, recorded only: calcified band " + ", ".join(f"{v:+.1f}" for v in calc.phase_z)
      + "; paramagnetic cSS " + ", ".join(f"{v:+.1f}" for v in para.phase_z))
if not len(calc) or not (calc.phase_z < 0).all(): fail("calcified band not found or not diamagnetic")
if not len(para) or not (para.phase_z > 0).all(): fail("a paramagnetic cSS lesion read as diamagnetic")
c2, x2 = a2["candidates"], a2["excluded"]
kept = sorted(set(c2.lesion) - {0}); lost = sorted(set(c1.lesion) - {0} - set(kept) - {CALC})
if lost: fail(f"CSS_USE_MIMIC=1 lost paramagnetic cSS lesions {lost}")
e10 = x2[x2.reason.str.startswith("calcification")]
if not (e10.lesion == CALC).any() or (c2.lesion == CALC).any(): fail("calcified band not excluded as E10")
mbo = pd.concat([c1[c1.mb], a1["excluded"][a1["excluded"].mb]])
mb_note = (f"microbleed row detected as {len(mbo)} object(s), mb_cluster_n {list(mbo.mb_cluster_n)}"
           + (f", E11 with the flag on: {bool(x2[x2.mb].reason.str.startswith('cluster').any())}" if len(mbo) else
              " (not generated by the detector on this phantom - covered by the unit test)"))
print(f"CSS_USE_MIMIC=1: E10 excludes the calcified band ({e10.rule_version.iloc[0]}), paramagnetic cSS kept {kept}; {mb_note}")
print("RESULT mimic OK")
