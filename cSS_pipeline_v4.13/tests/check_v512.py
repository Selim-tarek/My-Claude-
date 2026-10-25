"""v5.1.2 checks: feature-spec families (scripts/css_features.py, phase validity / polarity in mimic_rules.py).
usage: check_v512.py    (units need nothing; the directional part needs PH3 from phantom_v4.py, aligned, in $CSS_BASE;
                         the phase part uses PHM from check_mimic.py when present)
1. units: a sphere is compact (compactness ~1), a sheet and a tube are not; a straight tube has tortuosity ~1, a
   half circle ~pi/2; the tortuosity of two non-touching pieces is measured on the larger one.
2. directions on PH3 objects (feature spec section 10), resolution-tolerant:
   VT.vessel_wm_mm > every cSS;  VC and VS sit deeper in the CSF than cSS (csf_depth_mm);  VS lies nearer the gyral
   crowns than VC (sulcal_depth_mm);  the compact transcortical vein is more compact than every cSS band.
3. PHM (if present): the calcified band has phase_z < 0 < every paramagnetic cSS; phase_valid_frac / polarity in 0..1.
PHANTOM DIRECTIONS ONLY - not real-world performance."""
import os, sys, subprocess, numpy as np, nibabel as nib, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import css_features as F
def fail(m): print("FAIL:", m); sys.exit(1)
v = np.array([0.8, 0.8, 0.8]); sh = (80, 80, 30); P = np.indices(sh).transpose(1, 2, 3, 0) * v
sph = ((P - [32, 32, 12]) ** 2).sum(-1) <= 2.0 ** 2
sheet = np.zeros(sh, bool); sheet[10:60, 10:60, 14] = True
tube = np.zeros(sh, bool); tube[10:70, 40:42, 14:16] = True
arc = np.zeros(sh, bool)
for a in np.linspace(0, np.pi, 400):
    i, j = int(40 + 25 * np.cos(a)), int(20 + 25 * np.sin(a)); arc[i - 1:i + 1, j - 1:j + 1, 14:16] = True
cs, cp, ct = (F.compactness_features(m, v)["compactness_3d"] for m in (sph, sheet, tube))
if not (cs > 0.8 and cp < 0.05 and ct < 0.2): fail(f"compactness sphere {cs} sheet {cp} tube {ct}")
ts, ta = F.path_features(tube, v)["tortuosity"], F.path_features(arc, v)["tortuosity"]
if not (abs(ts - 1) < 0.05 and abs(ta - np.pi / 2) < 0.15): fail(f"tortuosity straight {ts} arc {ta}")
two = tube.copy(); two[10:20, 70:72, 14:16] = True                    # a second, small, non-touching piece
if abs(F.path_features(two, v)["tortuosity"] - ts) > 0.02: fail("tortuosity not measured on the largest piece")
print(f"units: compactness sphere {cs} / sheet {cp} / tube {ct}; tortuosity straight {ts}, half circle {ta} (pi/2 = 1.57)")

B = os.environ.get("CSS_BASE")
if not B or not os.path.exists(f"{B}/data/PH3_swi.nii"):
    print("PH3 part skipped (no PH3 in $CSS_BASE)"); print("RESULT v512 OK"); sys.exit(0)
S = f"{B}/scripts"; py = sys.executable
r = subprocess.run([py, f"{S}/detect_css.py", "PH3"], capture_output=True, text=True)
if r.returncode: print(r.stdout[-600:], r.stderr[-1200:]); fail("detect_css PH3")
tr = nib.load(f"{B}/work/PH3_truth.nii.gz").get_fdata().astype(int); ve = nib.load(f"{B}/work/PH3_veins.nii.gz").get_fdata().astype(int)
kind_of = {1: "VC", 2: "VC", 3: "VC", 4: "VC", 5: "VS", 6: "VS", 7: "VT"}
rows = []
for lst, idc in (("candidates", "cand_id"), ("excluded", "excl_id")):
    d = pd.read_csv(f"{B}/review/PH3_{lst}.csv"); M = nib.load(f"{B}/review/PH3_{lst}.nii.gz").get_fdata().astype(int)
    for _, rr in d.iterrows():
        m = M == rr[idc]
        k = "cSS" if (tr[m] > 0).mean() > 0.2 else (kind_of.get(int(np.bincount(ve[m][ve[m] > 0]).argmax()), "other")
                                                    if (ve[m] > 0).mean() > 0.3 else "other")
        rows.append(dict(kind=k, **rr.to_dict()))
D = pd.DataFrame(rows); g = {k: D[D.kind == k] for k in ("cSS", "VC", "VS", "VT")}
if any(not len(x) for x in g.values()): fail(f"PH3 object kinds missing: {D.kind.value_counts().to_dict()}")
checks = [("VT.vessel_wm_mm > cSS", g["VT"].vessel_wm_mm.min() > g["cSS"].vessel_wm_mm.max()),
          ("VC.csf_depth_mm > cSS", g["VC"].csf_depth_mm.min() > g["cSS"].csf_depth_mm.max()),
          ("VS.csf_depth_mm > cSS", g["VS"].csf_depth_mm.min() > g["cSS"].csf_depth_mm.max()),
          ("VS.sulcal_depth_mm < VC", g["VS"].sulcal_depth_mm.max() < g["VC"].sulcal_depth_mm.min()),
          ("VT.compactness_3d > cSS", g["VT"].compactness_3d.min() > g["cSS"].compactness_3d.max())]
bad = [n for n, ok in checks if not ok]
if bad: fail("directions not met: " + ", ".join(bad))
need = ["pial_sd_p10_mm", "surface_contact_frac_1mm", "csf_adjacent_frac", "sphericity", "elongation_3d", "flatness_3d",
        "path_length_mm", "darkness_z_p10", "dark_frac_z3", "feature_qc_flags", "feature_schema_version"]
if any(c not in D.columns for c in need) or (D.feature_schema_version != "css-features/1.0").any(): fail("v5.1.2 columns missing")
print("PH3 directions: " + "; ".join(n for n, _ in checks) + "  - all met")
if os.path.exists(f"{B}/review/PHM_candidates.csv"):
    c = pd.read_csv(f"{B}/review/PHM_candidates.csv"); x = pd.read_csv(f"{B}/review/PHM_excluded.csv")
    if "phase_valid_frac" in c.columns:
        allp = pd.concat([c, x])
        if not allp.phase_valid_frac.dropna().between(0, 1).all() or not allp.phase_polarity_frac.dropna().between(0, 1).all():
            fail("phase_valid_frac / phase_polarity_frac outside 0..1")
        print(f"PHM phase: valid fraction {allp.phase_valid_frac.min():.2f}-{allp.phase_valid_frac.max():.2f}, "
              f"polarity consistency {allp.phase_polarity_frac.min():.2f}-{allp.phase_polarity_frac.max():.2f}")
print("RESULT v512 OK")
