#!/usr/bin/env python3
"""Synthetic stress test across several host scans.
usage: stress_test.py --detector detect_css.py --tag v3 [--hosts P001,P002,P003,P004,P005]
                      [--depths 0.6,0.45,0.3] [--seeds 1] [--z -2.5]
Inserts synthetic cSS into each host (T1 anatomy is used when the host has it), runs the detector,
and measures sensitivity, worst
rank, how many real non-lesion candidates (veins/artifacts) outrank the lesions, and AUC."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import argparse, os, re, shutil, subprocess, sys, numpy as np, pandas as pd, nibabel as nib
from scipy import ndimage as ndi
ap = argparse.ArgumentParser()
ap.add_argument("--detector", default="detect_css.py"); ap.add_argument("--tag", default="run")
ap.add_argument("--hosts", default="P001,P002,P003,P004,P005")
ap.add_argument("--depths", default="0.6,0.45,0.3"); ap.add_argument("--seeds", default="1")
ap.add_argument("--z", default="-2.5")
ap.add_argument("--veins", default="4", help="synthetic vein decoys per synthetic scan (v3 generator)")
ap.add_argument("--radial", default="2", help="transcortical / medullary vein decoys per synthetic scan")
ap.add_argument("--legacy", action="store_true", help="v2 synthetic lesions (flood-filled, thicker)")
a = ap.parse_args()
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project")); S = f"{base}/scripts"
py = sys.executable
def run(*args):
    r = subprocess.run([py, *args], capture_output=True, text=True)
    if r.returncode: print(r.stdout[-500:], r.stderr[-1500:]); sys.exit("step failed: " + " ".join(args))
    return r.stdout
rows = []; feats = []; decoys = {"shown": 0, "excluded": 0, "not detected": 0}; why = {}; lost = {}
for h in a.hosts.split(","):
    clean = run(f"{S}/{a.detector}", h, a.z, "85", "3")
    nclean = int(re.search(r"-> (\d+) candidates", clean).group(1))
    # host anatomy: T1-derived labels (prep_anat.sh) when present, else SynthSeg on the SWI
    t1 = f"{base}/work/{h}_t1seg_swispace.nii.gz"
    use_t1 = os.path.exists(t1) and os.path.getmtime(t1) >= os.path.getmtime(f"{base}/data/{h}_swi.nii")
    if not use_t1:
        shutil.copy(f"{base}/synthseg/{h}_seg.nii.gz", f"{base}/synthseg/{h}S_seg.nii.gz")
        if os.path.exists(f"{base}/work/{h}S_t1seg_swispace.nii.gz"): os.remove(f"{base}/work/{h}S_t1seg_swispace.nii.gz")
    print(f"{h}: anatomy {'T1' if use_t1 else 'SWI-only'}", flush=True)
    for d in a.depths.split(","):
        for s in a.seeds.split(","):
            run(f"{S}/make_synthetic.py", h, "8", d, s, *(["--legacy"] if a.legacy else ["--veins", a.veins, "--radial", a.radial]))
            if use_t1:   # copied AFTER the synthetic SWI is written, so it counts as up to date
                shutil.copy(t1, f"{base}/work/{h}S_t1seg_swispace.nii.gz")
            run(f"{S}/align_seg.py", h + "S")
            run(f"{S}/{a.detector}", h + "S", a.z, "85", "3")
            res = re.search(r"RESULT (.*)", run(f"{S}/eval_synthetic.py", h + "S")).group(1)
            kv = dict(x.split("=") for x in res.split())
            row = dict(host=h, depth=float(d), seed=int(s), clean_candidates=nclean,
                       **{k: float(v) for k, v in kv.items()})
            rows.append(row)
            cdf = pd.read_csv(f"{base}/review/{h}S_candidates.csv")
            if len(cdf):
                tr = nib.load(f"{base}/work/{h}S_truth.nii.gz").get_fdata() > 0
                cm = nib.load(f"{base}/review/{h}S_candidates.nii.gz").get_fdata().astype(int)
                les = set(np.unique(cm[ndi.binary_dilation(tr, iterations=1)])) - {0}
                vp = f"{base}/work/{h}S_veins.nii.gz"
                vv = set(np.unique(cm[nib.load(vp).get_fdata() > 0])) - {0} if os.path.exists(vp) else set()
                cdf["is_lesion"] = cdf.cand_id.isin(les).astype(int); cdf["host"] = h; cdf["depth"] = float(d)
                cdf["is_decoy_vein"] = (cdf.cand_id.isin(vv) & ~cdf.cand_id.isin(les)).astype(int)
                # missed synthetic lesions: removed by which exclusion rule (or never detected)?
                tl = nib.load(f"{base}/work/{h}S_truth.nii.gz").get_fdata().astype(int)
                xm = nib.load(f"{base}/review/{h}S_excluded.nii.gz").get_fdata().astype(int)
                xr = pd.read_csv(f"{base}/review/{h}S_excluded.csv").set_index("excl_id")["reason"]
                # v4.14: components removed by the generation gates (size / slices / elongation / no cortex)
                dp_ = f"{base}/work/{h}S_dropped.nii.gz"
                dm = nib.load(dp_).get_fdata().astype(int) if os.path.exists(dp_) else np.zeros_like(xm)
                dr_ = pd.read_csv(f"{base}/work/{h}S_dropped.csv").set_index("drop_id")["reason"] \
                    if os.path.exists(f"{base}/work/{h}S_dropped.csv") else pd.Series(dtype=str)
                for t_ in range(1, int(tl.max()) + 1):
                    m_ = tl == t_
                    if (cm[m_] > 0).mean() >= 0.3: continue
                    if (xm[m_] > 0).any():
                        r_ = str(xr.get(int(np.bincount(xm[m_][xm[m_] > 0]).argmax()), "?")).split(":")[0].split(" (")[0]
                    elif (dm[m_] > 0).any():
                        r_ = "dropped: " + str(dr_.get(int(np.bincount(dm[m_][dm[m_] > 0]).argmax()), "?"))
                    else:
                        r_ = "not detected"
                    lost[r_] = lost.get(r_, 0) + 1
                # what happened to each decoy vein: still shown for review, excluded (why), or never detected
                if os.path.exists(vp):
                    vm = nib.load(vp).get_fdata().astype(int)
                    for v_ in range(1, int(vm.max()) + 1):
                        m_ = vm == v_
                        if (cm[m_] > 0).any(): decoys["shown"] += 1
                        elif (xm[m_] > 0).any():
                            decoys["excluded"] += 1
                            r_ = str(xr.get(int(np.bincount(xm[m_][xm[m_] > 0]).argmax()), "?")).split(":")[0].split(" (")[0]
                            why[r_] = why.get(r_, 0) + 1
                        else: decoys["not detected"] += 1
                feats.append(cdf)
            print(f"{h} depth {d} seed {s}: found {int(row['found'])}/{int(row['total'])}, "
                  f"worst #{int(row['worst'])}, non-lesions above {int(row['above'])}, AUC {row['auc']:.2f}", flush=True)
df = pd.DataFrame(rows)
os.makedirs(f"{base}/results", exist_ok=True)
df.to_csv(f"{base}/results/stress_{a.tag}.csv", index=False)
print(f"\n===== SUMMARY ({a.tag}, detector {a.detector}) =====")
for d, g in df.groupby("depth", sort=False):
    print(f"darkness {d}: found {int(g.found.sum())}/{int(g.total.sum())} "
          f"({100*g.found.sum()/g.total.sum():.0f}%)  median worst rank #{g.worst.median():.0f}  "
          f"median non-lesions above {g.above.median():.0f}  mean AUC {g.auc.mean():.2f}")
print(f"clean-scan candidates per host: {df.groupby('host').clean_candidates.first().to_dict()}")
if lost: print(f"missed synthetic lesions: " + ", ".join(f"{v} {k}" for k, v in lost.items()))
if sum(decoys.values()):
    print(f"decoy veins: {decoys['shown']} still shown, {decoys['excluded']} excluded "
          f"({', '.join(f'{v} {k}' for k, v in why.items()) or '-'}), {decoys['not detected']} not detected")
print(f"saved: {base}/results/stress_{a.tag}.csv")
if feats:
    F = pd.concat(feats, ignore_index=True)
    F.to_csv(f"{base}/results/stress_{a.tag}_features.csv", index=False)
    def auc(x, y):   # P(feature of a lesion > feature of a non-lesion)
        x, y = np.asarray(x, float), np.asarray(y, float); x, y = x[~np.isnan(x)], y[~np.isnan(y)]
        if not len(x) or not len(y): return float("nan")
        return float(((x[:, None] > y[None, :]).mean() + 0.5 * (x[:, None] == y[None, :]).mean()))
    print("\nFeature check: does it separate synthetic cSS from the scans' real veins/artifacts?")
    print("(AUC 0.5 = useless; >0.7 or <0.3 = useful; <0.5 means lower values indicate cSS)")
    for f in ["darkness_z", "volume_mm3", "cortex_frac", "cortex_dist_mm", "branch_per10mm",
              "surface_gradient", "surface_contact", "elongation", "pial_dist_mm", "bank_frac",
              "tube_ratio", "surface_alignment", "tram_frac", "vein_tree_mm", "mirror_dark_frac",
              "vessel_ext_mm", "vessel_wm_mm", "vessel_run_mm", "axis_normal", "score_v3", "score_v4", "score"]:
        if f in F.columns:
            L, N = F[F.is_lesion == 1][f], F[F.is_lesion == 0][f]
            V = F[F.get("is_decoy_vein", 0) == 1][f] if "is_decoy_vein" in F else N.iloc[0:0]
            print(f"  {f:18s} cSS median {L.median():7.2f} | others median {N.median():7.2f} | AUC {auc(L, N):.2f}"
                  + (f" | vs decoy veins {auc(L, V):.2f} (n={len(V)})" if len(V) else ""))
