"""v4.18 checks: ridge scales in mm, mirror lookup restricted to the brain, skeleton length with diagonal steps.
usage: check_v418.py    (needs PH1 and PH3 detected in $CSS_BASE)"""
import os, sys, json, numpy as np, pandas as pd
B = os.environ["CSS_BASE"]; sys.path.insert(0, f"{B}/scripts")
from css_common import skeleton_length
def fail(m): print("FAIL:", m); sys.exit(1)

# 1. skeleton length: straight, diagonal, L-corner (no double count at the corner)
z = np.zeros((20, 20), bool)
a = z.copy(); a[5, 2:12] = 1
b = z.copy(); b[np.arange(10), np.arange(10)] = 1
c = z.copy(); c[5, 2:8] = 1; c[5:10, 7] = 1
got = [skeleton_length(a, (0.5, 0.5)), skeleton_length(b, (0.5, 0.5)), skeleton_length(c, (1, 1))]
want = [4.5, 9 * np.hypot(0.5, 0.5), 9.0]
if not np.allclose(got, want): fail(f"skeleton length {got} != {want}")
print(f"   skeleton length: straight {got[0]:.2f}, diagonal {got[1]:.2f} (old rule {10 * 0.5:.2f}), corner {got[2]:.1f} mm")

# 2. ridge scales: the same mm on every protocol (PH1 0.5 mm, PH3 0.8 mm in-plane)
for s, inpl in (("PH1", 0.5), ("PH3", 0.8)):
    p = json.load(open(f"{B}/review/{s}_manifest.json"))["detect"]["parameters"]
    if not np.allclose(np.array(p["sato_sigmas_px"]) * inpl, p["sato_sigmas_mm"], atol=0.01):
        fail(f"{s}: sato px {p['sato_sigmas_px']} x {inpl} mm != {p['sato_sigmas_mm']} mm")
    print(f"   {s}: ridge sigma {p['sato_sigmas_mm']} mm = {p['sato_sigmas_px']} px")

# 3. mirror: a fraction in [0, 1] or NaN (most mirror points outside the brain), never forced dark by background
m = pd.concat([pd.read_csv(f"{B}/review/PH3_{k}.csv") for k in ("candidates", "excluded")]).mirror_dark_frac
if ((m < 0) | (m > 1)).any(): fail("mirror_dark_frac outside [0, 1]")
print(f"   PH3 mirror_dark_frac: {int(m.isna().sum())} NaN, median {m.median():.2f}")
# 4. end types on synthetic geometry: a 10 mm candidate whose far end stops / runs on / forks / reaches the midline
import css_vessel as cv, nibabel as nib
sh = (80, 60, 8); vx = (0.5, 0.5, 0.5)
for kind in ("stops", "runs_on", "branches", "midline"):
    cm = np.zeros(sh, bool); ves = np.zeros(sh, bool); mid = np.zeros(sh, bool)
    cm[10:30, 30, 4] = True
    if kind != "stops": ves[31:60, 30, 4] = True
    if kind == "branches":
        for k in range(16): ves[45 + k, 30 + k, 4] = True; ves[45 + k, 30 - k, 4] = True
    if kind == "midline": mid[55:] = True
    o = cv.local_tracking(cm, ves, vx, 20.0, dark=np.zeros(sh, bool), midline=mid)
    if (o["end_type_max"], o["end_type_min"]) != (kind, "stops"): fail(f"end type {kind}: {o['end_type_max']}/{o['end_type_min']}")
print("   end types: stops / runs_on / branches / midline recovered on synthetic geometry")

# 5. review tier on PH3 (phantom only - not a validation): every true lesion in tier 1
c = pd.read_csv(f"{B}/review/PH3_candidates.csv")
C = nib.load(f"{B}/review/PH3_candidates.nii.gz").get_fdata().astype(int)
t = nib.load(f"{B}/work/PH3_truth.nii.gz").get_fdata() > 0
c["truth"] = [t[C == i].mean() > 0.2 for i in c.cand_id]
if (c[c.truth].review_tier != 1).any(): fail(f"PH3 true lesion in tier 2: {c[c.truth & (c.review_tier != 1)].cand_id.tolist()}")
print(f"   PH3 tiers: lesions {c[c.truth].review_tier.tolist()}, others {c[~c.truth].review_tier.tolist()}")
# 6. radius (half-width) along the centre line: 1 / 3 / 5-voxel-wide bands at 0.5 mm -> 0.25 / 0.75 / 1.25 mm
from css_common import radius_profile
got = []
for w in (1, 3, 5):
    m = np.zeros((40, 40, 3), bool); m[10:30, 20 - w // 2: 21 + w // 2, 1] = True
    got.append(float(np.median(radius_profile(m, (0.5, 0.5, 1.0)))))
if not np.allclose(got, [0.25, 0.75, 1.25]): fail(f"radius profile {got}")
cols = ["cand_radius_med_mm", "cand_radius_p95_mm", "cand_radius_cv", "net_radius_med_mm", "net_radius_p95_mm"]
if not set(cols) <= set(c.columns): fail("radius columns missing")
print(f"   radius: bands {got} mm; PH3 lesions median {c[c.truth].cand_radius_med_mm.median():.2f} mm, "
      f"others {c[~c.truth].cand_radius_med_mm.median():.2f} mm (phantom only)")
print("RESULT v4.18 OK")
