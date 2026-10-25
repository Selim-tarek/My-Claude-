#!/usr/bin/env python3
"""How stable is the detector when the SAME scan looks a little different? (v4.18, idea from BrainIAC's
perturbation test). Runs detect_css.py on perturbed copies of one subject and compares each run with an
unperturbed copy. Changes nothing in the subject's own outputs; reported only.
usage: robustness.py SUBJECT [--perturb copy,inplane:0.86,slice:2.0,bias:0.3,noise:0.05,gamma:0.8,gamma:1.25,shift:0]
                             [--old OLD_candidates.nii.gz --ids 3,4,5] [--match-mm 5] [--top 10] [--keep]
Perturbations (applied inside the brain only, background stays 0):
  copy           unperturbed copy: the reference every other run is compared with (same inputs: no FLAIR,
                 no drawn ICH - copies never use them)
  inplane:MM     resample in-plane to MM voxels (anti-aliased when coarser) - another protocol
  slice:MM       resample the slice axis to MM slices
  bias:C         smooth multiplicative bias field exp(C * f), f in [-1, 1] across the brain
  noise:F        Gaussian noise, SD = F x the median brain intensity
  gamma:G        contrast change: I_max (I / I_max)^G
  shift:AXIS     image moved 1 voxel along array AXIS against unchanged labels (1-voxel misregistration)
Per run: candidates, excluded, reference candidates found again (centroid within --match-mm, one-to-one),
new candidates, Spearman correlation of scores of matched candidates, top-k overlap, and - with --old/--ids,
e.g. P006's v3 map - how many known lesions are still shown (as check_known.py).
Writes review/ID_robustness.csv. Copies are named ID_R<k>; their files are removed unless --keep.
Each run is a full detection: on a real scan expect several minutes per perturbation."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import argparse, glob, os, re, shutil, subprocess, sys, numpy as np, pandas as pd, nibabel as nib
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("subject")
ap.add_argument("--perturb", default="copy,inplane:0.86,slice:2.0,bias:0.3,noise:0.05,gamma:0.8,gamma:1.25,shift:0")
ap.add_argument("--old", default=""); ap.add_argument("--ids", default="")
ap.add_argument("--match-mm", type=float, default=5.0); ap.add_argument("--top", type=int, default=10)
ap.add_argument("--keep", action="store_true"); ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args(); subj = a.subject; base = _base(); S = os.path.dirname(os.path.abspath(__file__))
py = sys.executable

def run(*args):
    r = subprocess.run([py, *args], capture_output=True, text=True)
    if r.returncode: print(r.stdout[-600:], r.stderr[-1500:]); sys.exit("step failed: " + " ".join(args))
    return r.stdout

swi = nib.load(f"{base}/data/{subj}_swi.nii"); I0 = swi.get_fdata().astype(np.float32); A0 = swi.affine
vox0 = np.array(swi.header.get_zooms()[:3], float)
brain0 = I0 > 0
t1 = f"{base}/work/{subj}_t1seg_swispace.nii.gz"
lab_src = t1 if os.path.exists(t1) and os.path.getmtime(t1) >= os.path.getmtime(f"{base}/data/{subj}_swi.nii") \
    else f"{base}/synthseg/{subj}_seg.nii.gz"
use_t1 = lab_src == t1
rng = np.random.default_rng(a.seed)

def regrid(I, brain, new_vox):
    """resample to new voxel sizes on the same field of view (linear image, nearest mask); Gaussian
    anti-aliasing along axes that get coarser (FWHM^2 difference, FWHM = voxel size)"""
    new_vox = np.asarray(new_vox, float); f = new_vox / vox0
    sig = np.sqrt(np.clip(new_vox ** 2 - vox0 ** 2, 0, None)) / 2.355 / vox0
    Is = ndi.gaussian_filter(I, sig) if (sig > 0).any() else I
    shape = np.maximum(1, np.round(np.array(I.shape) / f)).astype(int)
    # new voxel j has its centre at old index (j + 0.5) * f - 0.5 (same outer edges)
    off = 0.5 * f - 0.5
    In = ndi.affine_transform(Is, np.diag(f), offset=off, output_shape=tuple(shape), order=1, mode="nearest")
    Bn = ndi.affine_transform(brain.astype(np.float32), np.diag(f), offset=off, output_shape=tuple(shape),
                              order=0, mode="constant") > 0.5
    An = A0.copy(); An[:3, :3] = A0[:3, :3] * f[None, :]; An[:3, 3] = A0[:3, :3] @ off + A0[:3, 3]
    return np.where(Bn, np.clip(In, 1e-3, None), 0).astype(np.float32), An

def perturb(kind, val):
    I, A = I0.copy(), A0.copy(); b = brain0
    med = float(np.median(I0[b]))
    if kind == "copy":
        pass
    elif kind == "inplane":
        I, A = regrid(I0, b, [val, val, vox0[2]])
    elif kind == "slice":
        I, A = regrid(I0, b, [vox0[0], vox0[1], val])
    elif kind == "bias":
        c = np.argwhere(b); lo, hi = c.min(0), c.max(0)
        g = [np.linspace(-1, 1, n) for n in I0.shape]
        w = rng.normal(size=3); w /= np.linalg.norm(w)
        # smooth field: random direction ramp + a gentle quadratic, scaled to [-1, 1] inside the brain
        X = np.meshgrid(*[(gi - gi[(l + h) // 2]) for gi, l, h in zip(g, lo, hi)], indexing="ij")
        fld = w[0] * X[0] + w[1] * X[1] + w[2] * X[2] + 0.5 * (X[0] ** 2 + X[1] ** 2 - X[2] ** 2)
        fb = fld[b]; fld = 2 * (fld - fb.min()) / (fb.max() - fb.min() + 1e-9) - 1
        I = np.where(b, I0 * np.exp(val * fld), 0)
    elif kind == "noise":
        I = np.where(b, np.clip(I0 + rng.normal(0, val * med, I0.shape), 1e-3, None), 0)
    elif kind == "gamma":
        mx = float(I0[b].max()); I = np.where(b, mx * (np.clip(I0, 0, None) / mx) ** val, 0)
    elif kind == "shift":
        ax = int(val); I = np.roll(I0, 1, axis=ax)
        sl = [slice(None)] * 3; sl[ax] = 0; I[tuple(sl)] = 0
    else:
        sys.exit(f"unknown perturbation {kind}")
    return I.astype(np.float32), A

def world(df, A):
    if not len(df): return np.zeros((0, 3))
    return (A @ np.c_[df[["vox_i", "vox_j", "vox_k"]].values, np.ones(len(df))].T)[:3].T

def known_kept(name):
    if not (a.old and a.ids): return np.nan
    old = nib.load(os.path.expanduser(a.old))
    cand = nib.load(f"{base}/review/{name}_candidates.nii.gz")
    from nibabel.processing import resample_from_to
    O = resample_from_to(old, cand, order=0).get_fdata().astype(int)
    C = cand.get_fdata().astype(int); k = 0
    for i in [int(x) for x in a.ids.split(",") if x.strip()]:
        m = ndi.binary_dilation(O == i, iterations=2)
        k += int(bool(set(np.unique(C[m])) - {0}))
    return k

def cleanup(name):
    for p in glob.glob(f"{base}/*/{name}_*"):
        if os.path.isfile(p): os.remove(p)

items = [(p.split(":")[0], float(p.split(":")[1]) if ":" in p else 0.0) for p in a.perturb.split(",") if p]
if items[0][0] != "copy": items = [("copy", 0.0)] + items
print(f"{subj}: robustness - {len(items)} detection runs on perturbed copies (anatomy {'T1' if use_t1 else 'SWI-only'}; "
      f"copies never use FLAIR or a drawn ICH mask)", flush=True)
rows, ref = [], None
for k, (kind, val) in enumerate(items):
    name = f"{subj}_R{k}"; label = kind if kind == "copy" else f"{kind}:{val:g}"
    I, A = perturb(kind, val)
    nib.save(nib.Nifti1Image(I, A), f"{base}/data/{name}_swi.nii")
    os.makedirs(f"{base}/work", exist_ok=True); os.makedirs(f"{base}/synthseg", exist_ok=True)
    shutil.copy(lab_src, f"{base}/work/{name}_t1seg_swispace.nii.gz" if use_t1 else f"{base}/synthseg/{name}_seg.nii.gz")
    run(f"{S}/align_seg.py", name)
    st = subprocess.run([py, f"{S}/detect_css.py", name], capture_output=True, text=True)
    if st.returncode:
        print(f"  {label:12s} detection failed / stopped: {st.stdout.strip().splitlines()[-1] if st.stdout.strip() else st.stderr[-300:]}")
        rows.append(dict(perturbation=label, status="failed")); continue
    c = pd.read_csv(f"{base}/review/{name}_candidates.csv"); x = pd.read_csv(f"{base}/review/{name}_excluded.csv")
    w = world(c, A)
    r = dict(perturbation=label, status="ok", voxel_mm=" x ".join(f"{v:.2f}" for v in nib.load(f"{base}/data/{name}_swi.nii").header.get_zooms()[:3]),
             candidates=len(c), excluded=len(x), known_kept=known_kept(name))
    if ref is None:
        ref = (c, w); r.update(found_again=len(c), new=0, spearman=1.0, top_overlap=min(a.top, len(c)))
    else:
        rc, rw = ref; pairs = []
        if len(rc) and len(c):
            d, j = cKDTree(w).query(rw)
            used = set()
            for i_ in np.argsort(d):                       # one-to-one, closest pairs first
                if d[i_] <= a.match_mm and j[i_] not in used:
                    pairs.append((i_, j[i_])); used.add(j[i_])
        r["found_again"] = len(pairs); r["new"] = len(c) - len(pairs)
        if len(pairs) >= 3:
            s1 = rc.score.values[[p[0] for p in pairs]]; s2 = c.score.values[[p[1] for p in pairs]]
            r["spearman"] = round(float(spearmanr(s1, s2).correlation), 2)
        else:
            r["spearman"] = np.nan
        topr = set(range(min(a.top, len(rc)))); topc = set(range(min(a.top, len(c))))
        r["top_overlap"] = sum(1 for p in pairs if p[0] in topr and p[1] in topc)
    rows.append(r)
    print(f"  {label:12s} {r['voxel_mm']:>18s} mm  candidates {r['candidates']:4d}  excluded {r['excluded']:4d}  "
          f"found again {r['found_again']:4d}/{len(ref[0])}  new {r['new']:4d}  spearman {r['spearman']}  "
          f"top-{a.top} {r['top_overlap']}" + ("" if np.isnan(r["known_kept"]) else
                                             f"  known kept {r['known_kept']}/{len(a.ids.split(','))}"), flush=True)
    if not a.keep: cleanup(name)
out = pd.DataFrame(rows); out.insert(0, "subject", subj)
out.to_csv(f"{base}/review/{subj}_robustness.csv", index=False)
print(f"{subj}: -> review/{subj}_robustness.csv (reported only; the subject's own outputs are unchanged)")
