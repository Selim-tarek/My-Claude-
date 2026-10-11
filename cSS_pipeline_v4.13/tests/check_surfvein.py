"""v4.16 checks on PH6 (PH3 + a surface-following, branching vein draining to a sagittal sinus):
local vein tracking (detect_css.py / css_vessel.py) and pial-surface features (surface_features.py).
usage: check_surfvein.py    (needs PH6 from phantom_surfvein.py, aligned, in $CSS_BASE)
PHANTOM REGRESSION CHECKS ONLY: they show that the features measure what they are meant to on constructed
anatomy. They do not validate any threshold for real data."""
import os, sys, subprocess, numpy as np, nibabel as nib, pandas as pd
B = os.environ["CSS_BASE"]; S = f"{B}/scripts"; py = sys.executable
def fail(m): print("FAIL:", m); sys.exit(1)
def run(*a, ok=True):
    r = subprocess.run([py, *a], capture_output=True, text=True)
    if ok and r.returncode: print(r.stdout[-800:], r.stderr[-1500:]); fail("step failed: " + " ".join(a))
    return r
run(f"{S}/detect_css.py", "PH6")
tr = nib.load(f"{B}/work/PH6_truth.nii.gz").get_fdata().astype(int)
sv = nib.load(f"{B}/work/PH6_surfvein.nii.gz").get_fdata().astype(int)
objs = []
for lst, idc in (("candidates", "cand_id"), ("excluded", "excl_id")):
    d = pd.read_csv(f"{B}/review/PH6_{lst}.csv"); M = nib.load(f"{B}/review/PH6_{lst}.nii.gz").get_fdata().astype(int)
    for _, r in d.iterrows():
        m = M == r[idc]
        kind = "cSS" if (tr[m] > 0).mean() > 0.2 else (
            {1: "surface vein", 2: "surface vein", 3: "sinus"}[int(np.bincount(sv[m][sv[m] > 0]).argmax())]
            if (sv[m] > 0).mean() > 0.5 else "other")
        objs.append(dict(list=lst[:-1] if lst == "candidates" else "excluded", obj_id=int(r[idc]), kind=kind, **r.to_dict()))
O = pd.DataFrame(objs)
found = [t for t in range(1, int(tr.max()) + 1)
         if (nib.load(f"{B}/review/PH6_candidates.nii.gz").get_fdata()[tr == t] > 0).mean() >= 0.3]
if len(found) < int(tr.max()) - 1: fail(f"cSS lost when the surface vein is added: {found}")
cols = ["kind", "list", "local_cont_max_mm", "terminal_ends", "local_branch_points", "score_v4"]
print(O[O.kind != "other"][cols].to_string(index=False))
c, v = O[O.kind == "cSS"], O[O.kind.isin(["surface vein", "sinus"])]
if not len(v): fail("surface vein / sinus not generated at all")
# 1. local tracking: the vein network runs on beyond the candidate (into the sinus / the other vein part)
vv = O[O.kind == "surface vein"]
if not len(vv) or not (vv.local_cont_max_mm >= 10).all() or not (vv.terminal_ends <= 1).all():
    fail("the surface vein should continue (into the sinus) beyond one end")
# sinus pieces (the 40 mm merge cap splits the sinus into several candidates) are reported, not asserted
if (c.local_cont_max_mm >= 10).any(): fail("a phantom cSS candidate continues >= 10 mm")
mv = O[(O.kind == "surface vein") & (O.cand_geodesic_mm == O[O.kind == "surface vein"].cand_geodesic_mm.max())]
if (mv.local_branch_points < 1).all(): fail("the branching surface vein shows no junction")
print(f"vein tracking: cSS continuation <= {c.local_cont_max_mm.max():.1f} mm, surface vein >= "
      f"{vv.local_cont_max_mm.min():.1f} mm (sinus pieces {', '.join(f'{x:.1f}' for x in O[O.kind == 'sinus'].local_cont_max_mm)}); "
      f"main vein junctions {int(mv.local_branch_points.max())}")

# 2. surface features: transform direction found from the labels, explicit transform gives the same values
run(f"{os.path.dirname(os.path.abspath(__file__))}/phantom_surfaces.py", "PH6", "6")
out = run(f"{S}/surface_features.py", "PH6").stdout
if "(LTA maps SWI->T1)" not in out: fail(f"wrong transform direction chosen: {out}")
F1 = pd.read_csv(f"{B}/review/PH6_surface_features.csv")
run(f"{S}/surface_features.py", "PH6", "--surface-xfm", f"{B}/work/PH6_tkras_to_swi.txt")
F2 = pd.read_csv(f"{B}/review/PH6_surface_features.csv")
if not np.allclose(F1.select_dtypes("number").fillna(-9).values, F2.select_dtypes("number").fillna(-9).values, atol=0.02):
    fail("LTA route and explicit-transform route disagree")
# real FreeSurfer LTA format (type 0 VOX2VOX, volume info with "cras", as bbregister writes it) -> same features
run(f"{os.path.dirname(os.path.abspath(__file__))}/phantom_surfaces.py", "PH6", "6", "--lta-vox2vox")
out0 = run(f"{S}/surface_features.py", "PH6").stdout
F0 = pd.read_csv(f"{B}/review/PH6_surface_features.csv")
if "(LTA maps SWI->T1)" not in out0 or not np.allclose(F1.select_dtypes("number").fillna(-9).values,
                                                       F0.select_dtypes("number").fillna(-9).values, atol=0.02):
    fail(f"type-0 (vox2vox) LTA gives different results: {out0}")
bad = np.loadtxt(f"{B}/work/PH6_tkras_to_swi.txt"); bad[:3, 3] += 8.0
np.savetxt(f"{B}/work/PH6_bad_xfm.txt", bad)
if run(f"{S}/surface_features.py", "PH6", "--surface-xfm", f"{B}/work/PH6_bad_xfm.txt", ok=False).returncode == 0:
    fail("an 8 mm misregistered surface was accepted")
if len(F1) != len(O) or F1.surf_valid.min() < 1: fail("not every object projected onto the surface")
F = F1.merge(O[["list", "obj_id", "kind"]], on=["list", "obj_id"])
fc, fv = F[F.kind == "cSS"], F[F.kind.isin(["surface vein", "sinus"])]
if not (fc.surf_normal_dist_mm < 0).all() or not (fv.surf_normal_dist_mm > 0).all():
    fail("surface normal distance: cSS should be subpial (<0), the surface vein outside (>0)")
print(f"surface features: direction from labels OK, RAS2RAS / VOX2VOX LTA / explicit transform identical, 8 mm shift refused; "
      f"normal distance cSS {fc.surf_normal_dist_mm.max():.2f} mm max, surface vein {fv.surf_normal_dist_mm.min():.2f} mm min")
print(F[["kind", "list", "surf_normal_dist_mm", "surf_pieces", "surf_band_continuity", "geo_length_mm",
         "surf_dark_cont_mm", "surf_exit_paths"]][F.kind != "other"].to_string(index=False))
print("RESULT surfvein OK")
