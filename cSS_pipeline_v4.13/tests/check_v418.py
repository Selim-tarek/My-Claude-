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
print("RESULT v4.18 OK")
