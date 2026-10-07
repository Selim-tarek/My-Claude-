# All code — cSS pipeline v3.2

## scripts/run_css.sh

```bash
#!/bin/zsh
# One patient, start to finish (up to review).
# usage: run_css.sh SUBJECT_ID [path/to/swi.nii(.gz)] [--nostrip]
#   - with a path: imports the scan; skull-strips it with SynthStrip unless --nostrip
#     (use --nostrip only for scans that are already skull-stripped)
#   - set NOVIEW=1 to skip opening freeview (batch use)
set -e
S=$1; SRC=$2; STRIP=1
[[ "$3" == "--nostrip" ]] && STRIP=0
B=${CSS_BASE:-$HOME/css_project}; SC=$B/scripts
if [ -z "$S" ]; then
  echo "usage: run_css.sh SUBJECT_ID [path/to/swi.nii or .nii.gz] [--nostrip]"; exit 1
fi
mkdir -p $B/data $B/synthseg $B/work $B/review
if [ -n "$SRC" ]; then
  if [ $STRIP -eq 1 ]; then
    echo "[0/4] importing + skull-stripping (SynthStrip, CSF kept) $SRC"
    mri_convert "$SRC" $B/work/${S}_raw.nii.gz > /dev/null
    mri_synthstrip -i $B/work/${S}_raw.nii.gz -o $B/data/${S}_swi.nii \
      -m $B/work/${S}_brainmask.nii.gz > /dev/null
  else
    echo "[0/4] importing (no skull-strip) $SRC"
    mri_convert "$SRC" $B/data/${S}_swi.nii > /dev/null
  fi
fi
if [ ! -f $B/data/${S}_swi.nii ]; then
  echo "ERROR: $B/data/${S}_swi.nii not found - give the scan path as 2nd argument"; exit 1
fi
echo "[1/4] SynthSeg segmentation"
if [ ! -f $B/synthseg/${S}_seg.nii.gz ]; then
  mri_synthseg --i $B/data/${S}_swi.nii --o $B/synthseg/${S}_seg.nii.gz \
    --parc --robust --threads 8 > /dev/null
else
  echo "      (already done, skipping)"
fi
echo "[2/4] aligning segmentation to SWI"
python $SC/align_seg.py $S > /dev/null
echo "[3/4] detecting cSS candidates"
python $SC/detect_css.py $S -2.5 85 3
if [ -z "$NOVIEW" ]; then
  echo "[4/4] opening viewer"
  freeview -v $B/data/${S}_swi.nii \
    $B/review/${S}_candidates.nii.gz:colormap=lut:opacity=0.6 > /dev/null 2>&1 &
fi
echo ""
echo "NEXT: python $SC/export_review.py $S 40 | pbcopy    (paste into the workbook)"
echo "      review in freeview, then: python $SC/mark_css.py $S <cSS ids> [--ich <ids>]"
```

## scripts/run_all.sh

```bash
#!/bin/zsh
# Batch: run several patients without opening the viewer.
# usage: run_all.sh P006 P007 P008        (expects data/<ID>_swi.nii already imported)
B=${CSS_BASE:-$HOME/css_project}
for s in "$@"; do
  echo "======== $s"
  NOVIEW=1 $B/scripts/run_css.sh $s | head -1
done
```

## install.sh

```bash
#!/bin/zsh
# usage (from the package root):  zsh install.sh
B=$HOME/css_project; SRC=$(cd "$(dirname "$0")" && pwd)/scripts
mkdir -p $B/scripts $B/data $B/synthseg $B/work $B/review $B/results
if ls $B/scripts/*.py >/dev/null 2>&1; then
  BK=$B/scripts/backup_$(date +%Y%m%d_%H%M); mkdir -p $BK; cp $B/scripts/*.py $B/scripts/*.sh $BK/ 2>/dev/null
  echo "old scripts backed up to $BK"
fi
cp $SRC/*.py $SRC/*.sh $B/scripts/ && chmod +x $B/scripts/*.sh $B/scripts/*.py
grep -q 'css_project/scripts' ~/.zshrc || echo 'export PATH="$HOME/css_project/scripts:$PATH"' >> ~/.zshrc
echo "installed to $B/scripts  - open a new Terminal window"
```

## scripts/align_seg.py

```python
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib
from nibabel.processing import resample_from_to
from scipy import ndimage as ndi
subj = sys.argv[1]
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
swi = nib.load(f"{base}/data/{subj}_swi.nii")
seg = nib.load(f"{base}/synthseg/{subj}_seg.nii.gz")
seg_swi = resample_from_to(seg, swi, order=0).get_fdata().astype(np.int32)
n_before = (seg_swi > 0).sum()
seg_swi[~ndi.binary_fill_holes(swi.get_fdata() > 0)] = 0
out = nib.Nifti1Image(seg_swi, swi.affine, swi.header); out.set_data_dtype(np.int32)
nib.save(out, f"{base}/work/{subj}_seg_swispace.nii.gz")
print(f"{subj}: labelled voxels before mask {n_before:,}, after mask {(seg_swi>0).sum():,}")
```

## scripts/detect_css.py

```python
#!/usr/bin/env python3
"""cSS candidate detector v3.
usage: detect_css.py SUBJECT [z_thr=-2.5] [ridge_pct=85] [rim_mm=3]

Changes from v2 (ideas from published tools):
  - SHIVA-CMB: clip bright outliers at the 99.5th percentile before normalising
  - van Harten 2023 / hysteresis: weak voxels (z<-1.5) are kept when connected to a
    strong seed (z<-2.5) -> fewer broken-up lesions, better medium-strength detection
  - microbleednet / Cao 2026: record features that separate cSS from veins
    (cortex fraction, distance from cortex, branching) for review and a future classifier;
    only low cortex fraction (<0.35, typical of mid-sulcus veins) is penalised for now
  - long structures spanning >60% of the slab (sinus / fissure veins) are down-ranked
  - FIX: volumes were under-reported ~2.7x (regionprops already returns mm3 with spacing)
  - MGH cSS histology: 'surface_gradient' = how much darker the pial surface is than
    1-2.5 mm deeper (recorded only; does not change ranking until validated)
  - saves a 'dark' map so scoring can grow accepted lesions to their full extent
"""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
from skimage.filters import sato
from skimage.measure import regionprops
from skimage.morphology import skeletonize

subj   = sys.argv[1]
z_thr  = float(sys.argv[2]) if len(sys.argv) > 2 else -2.5
r_pct  = float(sys.argv[3]) if len(sys.argv) > 3 else 85
RIM_MM = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
Z_LOW, R_PCT_LOW = -2.0, 80           # hysteresis (weak) thresholds (v3.2: tighter, stops leaking into normal cortex)
GROW_Z = -1.5                          # 'dark' map used for region growing at scoring
MIN_MM3, MIN_SLICES, MIN_ELONG = 25.0, 2, 2.5   # 25 mm3 true volume ~= old setting
ARTIFACT_ZONES = ("orbitofrontal", "entorhinal", "temporalpole", "parahippocampal",
                  "fusiform", "inferiortemporal")
LH = [2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 17, 18, 26, 28]
RH = [41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60]

base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
swi = nib.load(f"{base}/data/{subj}_swi.nii")
I   = swi.get_fdata().astype(np.float32)
seg = nib.load(f"{base}/work/{subj}_seg_swispace.nii.gz").get_fdata().astype(np.int32)
vox = tuple(float(v) for v in swi.header.get_zooms()[:3]); vmm3 = float(np.prod(vox))
s2d = np.ones((3, 3, 1), bool); s3 = np.ones((3, 3, 3), bool)

brain  = ndi.binary_fill_holes(I > 0)
cortex = (seg >= 1000) | (seg == 3) | (seg == 42)
it = max(1, int(round(2.0 / vox[0])))
edge  = brain & ~ndi.binary_erosion(brain, structure=s2d, iterations=it)
shell = (ndi.binary_dilation(cortex, structure=s2d, iterations=2) | edge) & brain
surface = (~brain) | (seg == 24)
near_surface = ndi.binary_dilation(surface, structure=s2d, iterations=it)

# 1. robust normalisation
hi = np.percentile(I[brain], 99.5)
Ic = np.minimum(I, hi)
c = Ic[cortex & (I > 0)]
mu = float(np.median(c)); sd = float(1.4826 * np.median(np.abs(c - mu)) + 1e-6)
# v3.2 LOCAL reference: compare each voxel with cortex within ~10 mm, so regionally dark but
# normal cortex (e.g. iron-rich motor cortex) is not called abnormal
_sig = tuple(10.0 / v for v in vox)
_cw = ndi.gaussian_filter(cortex.astype(np.float32), _sig)
_cm = ndi.gaussian_filter((Ic * cortex).astype(np.float32), _sig)
local_mu = np.where(_cw > 0.05, _cm / np.maximum(_cw, 1e-6), mu)
z = (Ic - local_mu) / sd

# 2. skull-strip edge removal (3D, mm)
dist_mm = ndi.distance_transform_edt(brain, sampling=vox)
rim = brain & (dist_mm < RIM_MM)
If = Ic.copy(); If[~brain] = mu

# 3. curvilinear dark structures, slice by slice
R = np.zeros_like(I)
for k in range(I.shape[2]):
    if I[:, :, k].max() > 0:
        R[:, :, k] = sato(If[:, :, k], sigmas=[1, 2, 3], black_ridges=True)
zone = shell & ~rim
r_hi = np.percentile(R[shell], r_pct); r_lo = np.percentile(R[shell], R_PCT_LOW)
strong = zone & (z < z_thr) & (R > r_hi)
weak   = zone & (z < Z_LOW) & (R > r_lo)          # strong is a subset of weak

# 4. hemisphere map; hysteresis labelling separately in each hemisphere
Lm = np.isin(seg, LH) | ((seg >= 1000) & (seg < 2000))
Rm = np.isin(seg, RH) | (seg >= 2000)
hemiL = ndi.distance_transform_edt(~Lm, sampling=vox) <= ndi.distance_transform_edt(~Rm, sampling=vox)
dmid = np.where(hemiL, ndi.distance_transform_edt(hemiL, sampling=vox),
                ndi.distance_transform_edt(~hemiL, sampling=vox))
lab = np.zeros(I.shape, np.int32); n_raw = 0; nxt = 0
for side in (hemiL, ~hemiL):
    wl, nw = ndi.label(weak & side, structure=s3)
    n_raw += nw
    keep = np.unique(wl[strong & side]); keep = keep[keep > 0]
    remap = np.zeros(nw + 1, np.int32)
    remap[keep] = np.arange(1, len(keep) + 1) + nxt
    nz = wl > 0
    lab[nz] = remap[wl[nz]]
    nxt += len(keep)

os.makedirs(f"{base}/work", exist_ok=True)
nib.save(nib.Nifti1Image((zone & (z < GROW_Z)).astype(np.uint8), swi.affine),
         f"{base}/work/{subj}_dark.nii.gz")

# 5. features
cortex_dist = ndi.distance_transform_edt(~cortex, sampling=vox)      # 0 inside cortex
# depth below the pial surface for cortex voxels (distance to nearest non-cortex, non-WM voxel)
pial_depth = ndi.distance_transform_edt(cortex | np.isin(seg, [2, 41]), sampling=vox)
n_brain_sl = max(1, int(brain.any(axis=(0, 1)).sum()))
lut = {}
lutp = os.path.join(os.environ.get("FREESURFER_HOME", ""), "FreeSurferColorLUT.txt")
if os.path.exists(lutp):
    for line in open(lutp):
        t = line.split()
        if len(t) > 2 and t[0].isdigit():
            lut[int(t[0])] = t[1]

def branches_per_10mm(mm):
    """2D skeleton per slice: branch points per 10 mm of centre-line (veins branch, cSS rarely)"""
    bp = 0; length = 0.0
    k8 = np.ones((3, 3)); k8[1, 1] = 0
    for k in range(mm.shape[2]):
        sl = mm[:, :, k]
        if sl.sum() < 3: continue
        sk = skeletonize(sl)
        nb = ndi.convolve(sk.astype(int), k8, mode="constant")
        bp += int((sk & (nb >= 3)).sum()); length += sk.sum() * vox[0]
    return 10.0 * bp / max(length, 1.0)

objs = ndi.find_objects(lab)
rows = []
for p in regionprops(lab, spacing=vox):
    vol = float(p.area)                         # already mm3 because spacing is given
    if vol < MIN_MM3: continue
    sl = objs[p.label - 1]
    pad = tuple(slice(max(s.start - 3, 0), s.stop + 3) for s in sl)
    mm = lab[pad] == p.label
    n_sl = int(mm.any(axis=(0, 1)).sum())
    if n_sl < MIN_SLICES: continue
    elong = p.axis_major_length / max(p.axis_minor_length, 1e-3)
    if elong < MIN_ELONG: continue
    nb = seg[pad][ndi.binary_dilation(mm, iterations=2)]
    nb = nb[nb >= 1000]
    if nb.size == 0: continue
    own = seg[pad][mm]                           # hemisphere from own voxels
    n_l = int((np.isin(own, LH) | ((own >= 1000) & (own < 2000))).sum())
    n_r = int((np.isin(own, RH) | (own >= 2000)).sum())
    if n_l != n_r:
        side = nb[nb < 2000] if n_l > n_r else nb[nb >= 2000]
        if side.size: nb = side
    code = int(np.bincount(nb).argmax())
    region = lut.get(code, str(code))
    contact   = float(near_surface[pad][mm].mean())
    depth     = float(-z[pad][mm].mean())
    cfrac     = float(cortex[pad][mm].mean())
    cdist     = float(np.median(cortex_dist[pad][mm]))
    bdens     = branches_per_10mm(mm)
    # depth profile (MGH histology: cSS iron is densest in the outermost cortex):
    # darker at the pial surface than 1-2.5 mm deeper -> positive gradient
    nbh = ndi.binary_dilation(mm, structure=s2d, iterations=4) & cortex[pad]
    dp = pial_depth[pad]; zz = z[pad]
    outer = nbh & (dp <= 1.0); inner = nbh & (dp > 1.0) & (dp <= 2.5)
    sgrad = float(zz[inner].mean() - zz[outer].mean()) if outer.sum() >= 5 and inner.sum() >= 5 else float('nan')
    artz      = any(a in region for a in ARTIFACT_ZONES)
    midz      = bool(np.median(dmid[pad][mm]) < 5.0)
    vein_like = bool(cfrac < 0.35)   # branching is recorded only: thick cSS bands also give spiky skeletons
    longs     = bool(n_sl * vox[2] > 40.0)    # v3.2: > 40 mm top-to-bottom (sinus / large veins)
    score = depth * (0.25 + contact) * np.log1p(vol)
    for flag, f in ((artz, 0.5), (midz, 0.5), (vein_like, 0.6), (longs, 0.5)):
        if flag: score *= f
    c0 = np.argwhere(mm).mean(0) + [s.start for s in pad]
    rows.append(dict(hemi="L" if code < 2000 else "R", region=region, label=code,
                     volume_mm3=round(vol, 1), n_slices=n_sl, elongation=round(elong, 1),
                     surface_contact=round(contact, 2), darkness_z=round(depth, 1),
                     artifact_zone=int(artz), midline_zone=int(midz), score=round(score, 2),
                     vox_i=int(round(c0[0])), vox_j=int(round(c0[1])), vox_k=int(round(c0[2])),
                     cortex_frac=round(cfrac, 2), cortex_dist_mm=round(cdist, 2),
                     branch_per10mm=round(bdens, 2), surface_gradient=round(sgrad, 2),
                     vein_like=int(vein_like),
                     long_structure=int(longs), _lab=p.label))

cols = ["cand_id", "hemi", "region", "label", "volume_mm3", "n_slices", "elongation",
        "surface_contact", "darkness_z", "artifact_zone", "midline_zone", "score",
        "vox_i", "vox_j", "vox_k", "cortex_frac", "cortex_dist_mm", "branch_per10mm", "surface_gradient",
        "vein_like", "long_structure", "accept"]
keep = np.zeros(I.shape, np.int16)
if rows:
    df = pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)
    df.insert(0, "cand_id", np.arange(1, len(df) + 1))
    remap = np.zeros(lab.max() + 1, np.int16)
    remap[df["_lab"].values] = df["cand_id"].values
    keep = remap[lab]
    df["accept"] = ""
    df = df[cols]
else:
    df = pd.DataFrame(columns=cols)

os.makedirs(f"{base}/review", exist_ok=True)
df.to_csv(f"{base}/review/{subj}_candidates.csv", index=False)
nib.save(nib.Nifti1Image(keep, swi.affine), f"{base}/review/{subj}_candidates.nii.gz")
print(f"{subj}: {n_raw} raw components -> {len(df)} candidates "
      f"({int(df.artifact_zone.astype(int).sum()) if len(df) else 0} artifact zone, "
      f"{int(df.midline_zone.astype(int).sum()) if len(df) else 0} midline, "
      f"{int(df.vein_like.astype(int).sum()) if len(df) else 0} vein-like)\n")
if len(df):
    print(df.head(15)[["cand_id", "hemi", "region", "volume_mm3", "n_slices", "darkness_z",
                       "cortex_frac", "branch_per10mm", "vein_like", "score",
                       "vox_i", "vox_j", "vox_k"]].to_string(index=False))
```

## scripts/review_css.py

```python
#!/usr/bin/env python3
"""Interactive cSS review: shows each suspected candidate, you agree or disagree.
usage: review_css.py SUBJECT [--top 20] [--redo]

Keys (or click the buttons):
  y = cSS (agree)      v = vein      o = normal cortex      a = artifact
  i = near ICH         u = unsure    b = back one            [ / ] = move slice down/up
  q = finish and score
Decisions are saved after every key press, so you can quit and resume later.
At the end, accepted candidates are scored automatically (score_css.py)."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, subprocess, numpy as np, nibabel as nib, pandas as pd
import matplotlib
if os.environ.get("CSS_REVIEW_TEST"): matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.patches import Rectangle

subj = sys.argv[1]
top = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else 20
redo = "--redo" in sys.argv
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
csv = f"{base}/review/{subj}_candidates.csv"
calls_csv = f"{base}/review/{subj}_calls.csv"

swi = nib.as_closest_canonical(nib.load(f"{base}/data/{subj}_swi.nii"))      # RAS for display
cnd = nib.as_closest_canonical(nib.load(f"{base}/review/{subj}_candidates.nii.gz"))
I = swi.get_fdata().astype(np.float32); C = cnd.get_fdata().astype(int)
vox = swi.header.get_zooms()[:3]
br = I > 0
lo, hi = np.percentile(I[br], [1, 99]) if br.any() else (I.min(), I.max())
df = pd.read_csv(csv)
ids = df.cand_id.astype(int).tolist()[:top]
calls = {}
if os.path.exists(calls_csv) and not redo:
    calls = dict(pd.read_csv(calls_csv).astype({"cand_id": int}).values.tolist())
KEYS = {"y": "cSS", "v": "Vein", "o": "Normal", "a": "Artifact", "i": "Near ICH", "u": "Unsure"}
COL = {"cSS": "#2e7d32", "Vein": "#1565c0", "Normal": "#6d6d6d", "Artifact": "#ef6c00",
       "Near ICH": "#8e24aa", "Unsure": "#c9a400"}
HW = int(round(30 / vox[0]))          # 60 mm zoom window

class Reviewer:
    def __init__(self):
        self.pos = next((n for n, c in enumerate(ids) if c not in calls), 0)
        self.shift = 0
        self.fig = plt.figure(figsize=(14, 8.2))
        gs = self.fig.add_gridspec(2, 4, left=0.02, right=0.98, top=0.88, bottom=0.14, wspace=0.05, hspace=0.12)
        self.ov = self.fig.add_subplot(gs[:, 0])
        self.ax = [[self.fig.add_subplot(gs[r, c + 1]) for c in range(3)] for r in range(2)]
        self.btns = []
        labels = [("cSS (y)", "y"), ("Vein (v)", "v"), ("Normal (o)", "o"), ("Artifact (a)", "a"),
                  ("Near ICH (i)", "i"), ("Unsure (u)", "u"), ("◀ Back (b)", "b"),
                  ("Slice ↓ [", "["), ("Slice ↑ ]", "]"), ("Finish (q)", "q")]
        w = 0.094
        for n, (lab, k) in enumerate(labels):
            bax = self.fig.add_axes([0.02 + n * (w + 0.003), 0.03, w, 0.06])
            b = Button(bax, lab, color=COL.get(KEYS.get(k, ""), "#e0e0e0") if k in KEYS else "#e0e0e0",
                       hovercolor="#ffffff")
            if k in KEYS: b.label.set_color("white"); b.label.set_fontweight("bold")
            b.on_clicked(lambda e, k=k: self.key(k)); self.btns.append(b)
        self.fig.canvas.mpl_connect("key_press_event", lambda e: self.key(e.key))
        self.show()

    def show(self):
        if self.pos >= len(ids): return self.finish()
        cid = ids[self.pos]; r = df[df.cand_id == cid].iloc[0]
        m = C == cid
        if not m.any():
            self.pos += 1; return self.show()
        pts = np.argwhere(m); ci, cj = pts[:, 0].mean(), pts[:, 1].mean()
        ks = np.bincount(pts[:, 2]); kc = int(ks.argmax()) + self.shift
        kc = int(np.clip(kc, 1, I.shape[2] - 2))
        i0, i1 = int(max(ci - HW, 0)), int(min(ci + HW, I.shape[0]))
        j0, j1 = int(max(cj - HW, 0)), int(min(cj + HW, I.shape[1]))
        # overview
        self.ov.clear(); self.ov.imshow(I[:, :, kc].T, cmap="gray", origin="lower", vmin=lo, vmax=hi)
        self.ov.add_patch(Rectangle((i0, j0), i1 - i0, j1 - j0, fill=False, ec="yellow", lw=1.5))
        self.ov.set_title("whole slice (yellow = zoom)", fontsize=9); self.ov.axis("off")
        self.ov.text(2, I.shape[1] - 6, "L", color="yellow", fontsize=11); self.ov.text(I.shape[0] - 12, I.shape[1] - 6, "R", color="yellow", fontsize=11)
        # zoomed: top row raw, bottom row with outline; slices kc-1, kc, kc+1
        for c, k in enumerate((kc - 1, kc, kc + 1)):
            for row in (0, 1):
                a = self.ax[row][c]; a.clear()
                a.imshow(I[i0:i1, j0:j1, k].T, cmap="gray", origin="lower", vmin=lo, vmax=hi)
                if row == 1 and m[i0:i1, j0:j1, k].any():
                    a.contour(m[i0:i1, j0:j1, k].T.astype(float), levels=[0.5], colors="red", linewidths=1.2)
                a.set_xticks([]); a.set_yticks([])
                if row == 0: a.set_title(f"slice {k}" + ("  (centre)" if k == kc else ""), fontsize=9)
        self.ax[0][0].set_ylabel("raw SWI", fontsize=9); self.ax[1][0].set_ylabel("suspect outlined", fontsize=9)
        flags = [n for n, f in (("artifact zone", "artifact_zone"), ("midline", "midline_zone"), ("vein-like", "vein_like")) if int(r.get(f, 0) or 0)]
        done = sum(1 for c in ids if c in calls)
        prev = calls.get(cid)
        self.fig.suptitle(f"{subj}   candidate #{cid}  ({self.pos + 1}/{len(ids)}, {done} decided)   —   "
                          f"{r.hemi} {r.region}   {r.volume_mm3} mm³, {int(r.n_slices)} slices, darkness z {r.darkness_z}"
                          + (f"   [{', '.join(flags)}]" if flags else "")
                          + (f"\nyour previous call: {prev}" if prev else "\nIs the red outline cSS?"),
                          fontsize=11, color=COL.get(prev, "black"))
        self.fig.canvas.draw_idle()

    def key(self, k):
        if k in KEYS:
            calls[ids[self.pos]] = KEYS[k]; self.save(); self.pos += 1; self.shift = 0; self.show()
        elif k == "b":
            self.pos = max(self.pos - 1, 0); self.shift = 0; self.show()
        elif k == "[": self.shift -= 1; self.show()
        elif k == "]": self.shift += 1; self.show()
        elif k == "q": self.finish()

    def save(self):
        pd.DataFrame(sorted(calls.items()), columns=["cand_id", "call"]).to_csv(calls_csv, index=False)

    def finish(self):
        self.save(); plt.close(self.fig)

R = Reviewer()
if os.environ.get("CSS_REVIEW_TEST"):
    R.fig.savefig(f"{base}/review/{subj}_review_preview.png", dpi=80)
    for n, c in enumerate(ids): calls[c] = ["cSS", "Vein", "Near ICH", "Normal"][n % 4]
    R.save()
else:
    plt.show()

# ---- write decisions and score
accepted = sorted(c for c, v in calls.items() if v == "cSS")
ich = sorted(c for c, v in calls.items() if v == "Near ICH")
unsure = sorted(c for c, v in calls.items() if v == "Unsure")
print(f"\n{subj}: reviewed {len(calls)} of top {len(ids)}  |  cSS {accepted}  |  near-ICH {ich}"
      + (f"  |  unsure {unsure} (counted as NOT cSS - re-check with: review_css.py {subj})" if unsure else ""))
cmd = [sys.executable, f"{base}/scripts/mark_css.py", subj, ",".join(map(str, accepted)) or "none"]
if ich: cmd += ["--ich", ",".join(map(str, ich))]
subprocess.run(cmd, check=False)
print(f"decisions saved: {calls_csv}")
```

## scripts/export_review.py

```python
#!/usr/bin/env python3
"""Print top-N candidates as tab-separated rows for pasting into the workbook
('Candidate Review', column A).   usage: export_review.py P006 40 | pbcopy"""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, pandas as pd
subj = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 40
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
df = pd.read_csv(f"{base}/review/{subj}_candidates.csv").head(n)
cols = ["cand_id", "hemi", "region", "label", "volume_mm3", "n_slices", "elongation",
        "surface_contact", "darkness_z", "artifact_zone", "midline_zone", "score",
        "vox_i", "vox_j", "vox_k"]
for _, r in df.iterrows():
    print(subj + "\t" + "\t".join(str(r[c]) for c in cols))
```

## scripts/mark_css.py

```python
#!/usr/bin/env python3
"""Mark reviewed candidates and score.
usage: mark_css.py SUBJECT 1,4,9 [--ich 3,7]      (use 'none' when no cSS)
--ich: candidates you judged to be siderosis/hemosiderin CONNECTED TO a lobar ICH.
       They are excluded from the multifocality score (standard convention) but reported."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, subprocess, pandas as pd
if len(sys.argv) < 3:
    sys.exit("usage: mark_css.py SUBJECT 1,4,9 [--ich 3,7]   (or: none)")
subj, arg = sys.argv[1], sys.argv[2]
ich_arg = sys.argv[sys.argv.index("--ich") + 1] if "--ich" in sys.argv else ""
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
csv = f"{base}/review/{subj}_candidates.csv"
df = pd.read_csv(csv)
parse = lambda a: [] if a.lower() in ("", "none") else [int(x) for x in a.split(",") if x.strip()]
ids, ich = parse(arg), parse(ich_arg)
valid = set(df.cand_id)
bad = [i for i in ids + ich if i not in valid]
if bad: sys.exit(f"candidate ids not in list: {bad} (valid 1-{df.cand_id.max()})")
both = sorted(set(ids) & set(ich))
if both: sys.exit(f"ids marked both cSS and near-ICH: {both}")
df["accept"] = df.cand_id.isin(ids).astype(int)
df["near_ich"] = df.cand_id.isin(ich).astype(int)
df.to_csv(csv, index=False)
print(f"{subj}: accepted {len(ids)} of {len(df)} candidates {ids if ids else ''}"
      + (f"; near-ICH (not scored) {ich}" if ich else ""))
subprocess.run([sys.executable, f"{base}/scripts/score_css.py", subj], check=True)
```

## scripts/score_css.py

```python
#!/usr/bin/env python3
"""Score accepted candidates.  usage: score_css.py SUBJECT [--truth]
v2: adds van Harten-style seeded region growing - accepted candidates grow into connected
dark voxels (z<-1.5, inside the search zone, max 5 mm from the candidate) to give the
full lesion extent and a continuous cSS volume (no 0-4 ceiling effect).
Candidates marked near_ich=1 (siderosis connected to a lobar ICH) are excluded and reported."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, json, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi

subj = sys.argv[1]
MERGE_MM, ADJ_MM, GROW_MM = 3.0, 10.0, 5.0
use_truth = "--truth" in sys.argv
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
csv = f"{base}/review/{subj}_candidates.csv"
df  = pd.read_csv(csv)
img = nib.load(f"{base}/review/{subj}_candidates.nii.gz")
lab = img.get_fdata().astype(int); vox = tuple(float(v) for v in img.header.get_zooms()[:3])
vmm3 = float(np.prod(vox))

if use_truth:
    truth = nib.load(f"{base}/work/{subj}_truth.nii.gz").get_fdata() > 0
    df["accept"] = [int(((lab == c) & truth).any()) for c in df.cand_id]
    df.to_csv(csv, index=False)

df["accept"] = pd.to_numeric(df["accept"], errors="coerce")
if len(df) and df["accept"].isna().any():
    sys.exit(f"Unreviewed candidates (fill accept with 1 or 0): {df.loc[df.accept.isna(), 'cand_id'].tolist()}")

def min_dist(a, b):
    idx = np.argwhere(a | b); pad = int(np.ceil(ADJ_MM / min(vox))) + 1
    lo = np.maximum(idx.min(0) - pad, 0); hi = idx.max(0) + pad + 1
    sl = tuple(slice(l, h) for l, h in zip(lo, hi))
    return float(ndi.distance_transform_edt(~a[sl], sampling=vox)[b[sl]].min())

def groups(ids, D, thr):
    parent = {i: i for i in ids}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i in ids:
        for j in ids:
            if i < j and D[(i, j)] <= thr: parent[find(i)] = find(j)
    out = {}
    for i in ids: out.setdefault(find(i), []).append(i)
    return list(out.values())

acc = df[df.accept == 1] if len(df) else df
result = {"subject": subj}; total = total_foci = 0
for h in ["L", "R"]:
    ids = acc[acc.hemi == h].cand_id.astype(int).tolist() if len(acc) else []
    masks = {i: lab == i for i in ids}
    D = {(i, j): min_dist(masks[i], masks[j]) for i in ids for j in ids if i < j}
    foci = groups(ids, D, MERGE_MM); clusters = groups(ids, D, ADJ_MM)
    fpc = [sum(1 for f in foci if f[0] in c) for c in clusters]
    score = 0 if not ids else (1 if len(clusters) == 1 and fpc[0] <= 3 else 2)
    result[f"{h}_score"], result[f"{h}_foci"], result[f"{h}_clusters"] = score, len(foci), len(clusters)
    total += score; total_foci += len(foci)

# seeded region growing for full extent / volume
grown_vol = 0.0
seeds = np.isin(lab, acc.cand_id.astype(int).tolist()) if len(acc) else np.zeros(lab.shape, bool)
dark_p = f"{base}/work/{subj}_dark.nii.gz"
if seeds.any() and os.path.exists(dark_p):
    dark = nib.load(dark_p).get_fdata() > 0
    near = ndi.distance_transform_edt(~seeds, sampling=vox) <= GROW_MM
    grown = ndi.binary_propagation(seeds, structure=ndi.generate_binary_structure(3, 1),
                                   mask=(dark & near) | seeds)
    grown_vol = float(grown.sum() * vmm3)
    nib.save(nib.Nifti1Image(grown.astype(np.uint8), img.affine), f"{base}/review/{subj}_css_mask.nii.gz")

result["multifocality_0_4"] = total
result["n_foci_total"] = total_foci
result["category"] = "absent" if total_foci == 0 else ("focal" if total_foci <= 3 else "disseminated")
result["candidate_volume_mm3"] = round(float(acc.volume_mm3.sum()) if len(acc) else 0.0, 1)
result["grown_volume_mm3"] = round(grown_vol, 1)
result["regions"] = sorted(acc.region.unique().tolist()) if len(acc) else []
ich = df[df.near_ich == 1] if "near_ich" in df.columns and len(df) else df.iloc[0:0]
result["near_ich_candidates"] = int(len(ich))
result["near_ich_volume_mm3"] = round(float(ich.volume_mm3.sum()) if len(ich) else 0.0, 1)

json.dump(result, open(f"{base}/review/{subj}_score.json", "w"), indent=2)
summ = f"{base}/review/css_scores.csv"
row = pd.DataFrame([{k: v for k, v in result.items() if k != "regions"}])
if os.path.exists(summ):
    old = pd.read_csv(summ); old = old[old.subject != subj]
    row = pd.concat([old, row], ignore_index=True)
row.to_csv(summ, index=False)

print(f"\n{subj}  cSS multifocality score: {total}/4   ({result['category']})")
print(f"  left : score {result['L_score']}  foci {result['L_foci']}  clusters {result['L_clusters']}")
print(f"  right: score {result['R_score']}  foci {result['R_foci']}  clusters {result['R_clusters']}")
print(f"  volume: candidates {result['candidate_volume_mm3']} mm3, grown (full extent) {result['grown_volume_mm3']} mm3")
print(f"  regions: {', '.join(result['regions'])}")
if result["near_ich_candidates"]:
    print(f"  near-ICH siderosis (NOT in score): {result['near_ich_candidates']} candidates, "
          f"{result['near_ich_volume_mm3']} mm3")
```

## scripts/feature_report.py

```python
#!/usr/bin/env python3
"""Compare candidate features between reader calls (cSS vs Vein/Normal/Artifact) across all
reviewed subjects (review/*_calls.csv written by review_css.py).  usage: feature_report.py"""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, glob, numpy as np, pandas as pd
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
rows = []
for f in glob.glob(f"{base}/review/*_calls.csv"):
    s = os.path.basename(f).replace("_calls.csv", "")
    c = pd.read_csv(f); d = pd.read_csv(f"{base}/review/{s}_candidates.csv")
    m = d.merge(c, on="cand_id"); m["subject"] = s; rows.append(m)
if not rows: sys.exit("no reviewed cases yet (run review_css.py first)")
D = pd.concat(rows, ignore_index=True)
D = D[D.call.isin(["cSS", "Vein", "Normal", "Artifact"])]
print(f"{D.subject.nunique()} subject(s), {len(D)} decided candidates: {D.call.value_counts().to_dict()}\n")
def auc(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float); x, y = x[~np.isnan(x)], y[~np.isnan(y)]
    if not len(x) or not len(y): return float("nan")
    return float((x[:, None] > y[None, :]).mean() + 0.5 * (x[:, None] == y[None, :]).mean())
feats = ["darkness_z", "volume_mm3", "n_slices", "elongation", "surface_contact", "cortex_frac",
         "cortex_dist_mm", "branch_per10mm", "surface_gradient", "score"]
print(f"{'feature':18s} {'cSS':>8s} {'Vein':>8s} {'Normal':>8s} {'Artifact':>9s}   AUC cSS vs rest")
for f in feats:
    if f not in D: continue
    med = D.groupby("call")[f].median()
    a = auc(D[D.call == "cSS"][f], D[D.call != "cSS"][f])
    tag = "  <- useful" if a >= 0.75 or a <= 0.25 else ""
    print(f"{f:18s} " + " ".join(f"{med.get(k, np.nan):8.2f}" for k in ["cSS", "Vein", "Normal"]) +
          f" {med.get('Artifact', np.nan):9.2f}   {a:5.2f}{tag}")
print("\nAUC near 0.5 = no help; >0.75 = higher in cSS; <0.25 = lower in cSS.")
```

## scripts/make_synthetic.py

```python
#!/usr/bin/env python3
"""Insert synthetic cSS into a (presumed cSS-negative) host scan -> {HOST}S.
usage: make_synthetic.py HOST [n_lesions=8] [depth=0.6] [seed=1]
v2 (ideas from Cao et al. 2026): per-lesion random darkness, thickness (1-2 voxels),
length 10-30 mm over 2-4 slices, random edge blur, and added noise inside lesions,
so lesions are less 'clean' than before."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib
from scipy import ndimage as ndi

host  = sys.argv[1]
N     = int(sys.argv[2])   if len(sys.argv) > 2 else 8
depth = float(sys.argv[3]) if len(sys.argv) > 3 else 0.6
seed  = int(sys.argv[4])   if len(sys.argv) > 4 else 1
rng = np.random.default_rng(seed)
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))

img = nib.load(f"{base}/data/{host}_swi.nii"); I = img.get_fdata().astype(np.float32)
seg = nib.load(f"{base}/work/{host}_seg_swispace.nii.gz").get_fdata().astype(int)
vox = np.array(img.header.get_zooms()[:3], float)
s2 = np.ones((3, 3, 1), bool)
brain = ndi.binary_fill_holes(I > 0)
dist = ndi.distance_transform_edt(brain, sampling=vox)
cortex = seg >= 1000
CONVEX = [3, 8, 11, 15, 22, 24, 27, 28, 29, 30, 31]
convex = cortex & np.isin(seg % 1000, CONVEX)
sulcal = brain & ~cortex & ~np.isin(seg, [2, 41])
boundary = convex & ndi.binary_dilation(sulcal, structure=s2) & (dist > 4.0)
c = I[cortex & (I > 0)]; sd_c = 1.4826 * np.median(np.abs(c - np.median(c)))

nz = I.shape[2]
pts_all = np.argwhere(boundary)
pts_all = pts_all[(pts_all[:, 2] >= 3) & (pts_all[:, 2] < nz - 5)]
if len(pts_all) == 0: sys.exit("no admissible cortical boundary found")

truth, seeds, tries = np.zeros(I.shape, np.int16), [], 0
dark = np.zeros(I.shape, np.float32)
while len(seeds) < N and tries < 8000:
    tries += 1
    p = pts_all[rng.integers(len(pts_all))]
    if any(np.linalg.norm((p - q) * vox) < 25 for q in seeds): continue
    L_mm, nsl = rng.uniform(10, 30), int(rng.integers(2, 5))
    lesion, cur, ok = np.zeros(I.shape, bool), p[:2], True
    for dk in range(nsl):
        k = p[2] + dk
        if k >= nz: ok = False; break
        b = boundary[:, :, k]; cand = np.argwhere(b)
        if len(cand) == 0: ok = False; break
        d = np.linalg.norm((cand - cur) * vox[:2], axis=1); j = d.argmin()
        if d[j] > 4: ok = False; break
        start = tuple(cand[j]); cur = cand[j]
        target, seen, queue = int(L_mm / vox[0]), {start}, [start]
        while queue and len(seen) < target:
            x, y = queue.pop(0)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    nx_, ny_ = x + dx, y + dy
                    if 0 <= nx_ < b.shape[0] and 0 <= ny_ < b.shape[1] and b[nx_, ny_] and (nx_, ny_) not in seen:
                        seen.add((nx_, ny_)); queue.append((nx_, ny_))
        if len(seen) < 0.6 * target: ok = False; break
        for x, y in seen: lesion[x, y, k] = True
    if not ok: continue
    lesion = ndi.binary_dilation(lesion, structure=s2, iterations=int(rng.integers(1, 3))) & brain
    if (truth[lesion] > 0).any(): continue
    seeds.append(p); truth[lesion] = len(seeds)
    d_i = float(np.clip(depth + rng.uniform(-0.1, 0.1), 0.05, 0.9))
    w = ndi.gaussian_filter(lesion.astype(np.float32), sigma=(rng.uniform(0.6, 1.2),) * 2 + (0,))
    dark = np.maximum(dark, np.clip(w * 1.5, 0, 1) * d_i)

Is = I * (1 - dark)
Is += (rng.normal(0, 0.35 * sd_c, I.shape) * (dark > 0.05)).astype(np.float32)
Is[brain] = np.maximum(Is[brain], 1.0)        # keep the brain mask unchanged
Is[~brain] = I[~brain]

out = nib.Nifti1Image(Is.astype(np.float32), img.affine, img.header); out.set_data_dtype(np.float32)
nib.save(out, f"{base}/data/{host}S_swi.nii")
nib.save(nib.Nifti1Image(truth, img.affine), f"{base}/work/{host}S_truth.nii.gz")

lut = {}
lutp = os.path.join(os.environ.get("FREESURFER_HOME", ""), "FreeSurferColorLUT.txt")
if os.path.exists(lutp):
    for line in open(lutp):
        t = line.split()
        if len(t) > 2 and t[0].isdigit(): lut[int(t[0])] = t[1]
print(f"{host}S: inserted {len(seeds)} synthetic cSS lesions (depth ~{depth}, seed {seed})")
for i in range(1, len(seeds) + 1):
    m = truth == i; s = seg[m][seg[m] >= 1000]
    lab = int(np.bincount(s).argmax()) if s.size else 0
    ks = np.unique(np.argwhere(m)[:, 2]); c0 = np.argwhere(m).mean(0).round().astype(int)
    print(f"  lesion {i}: {lut.get(lab, lab)!s:30s} slices {ks.min()}-{ks.max()}  vox {c0[0]} {c0[1]} {c0[2]}")
```

## scripts/eval_synthetic.py

```python
#!/usr/bin/env python3
"""Score detector output against synthetic truth.  usage: eval_synthetic.py P005S
Reports sensitivity, worst rank, how many NON-lesion candidates (real veins/artifacts
of the host scan) rank above the worst lesion, and a ranking AUC."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
subj = sys.argv[1]
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
truth = nib.load(f"{base}/work/{subj}_truth.nii.gz").get_fdata().astype(int)
cand  = nib.load(f"{base}/review/{subj}_candidates.nii.gz").get_fdata().astype(int)
n = int(truth.max()); best, lesion_c = [], set()
for i in range(1, n + 1):
    ids = sorted(set(np.unique(cand[ndi.binary_dilation(truth == i, iterations=1)])) - {0})
    if ids: best.append(ids[0]); lesion_c |= set(ids)
    print(f"  lesion {i}: {'FOUND  (candidate #' + str(ids[0]) + ')' if ids else 'MISSED'}")
all_c = set(np.unique(cand)) - {0}
non = sorted(all_c - lesion_c)
worst = max(best) if best else 0
above = sum(1 for c in non if c < worst)
pairs = [(b < c) for b in best for c in non]
auc = float(np.mean(pairs)) if pairs else (1.0 if best else 0.0)
print(f"\n  sensitivity: {len(best)}/{n} = {100*len(best)/max(n,1):.0f}%   worst rank #{worst}")
print(f"  non-lesion candidates ranked above the worst lesion: {above}   ranking AUC: {auc:.2f}")
print(f"RESULT found={len(best)} total={n} worst={worst} above={above} auc={auc:.3f} ncand={len(all_c)}")
```

## scripts/stress_test.py

```python
#!/usr/bin/env python3
"""Synthetic stress test across several host scans.
usage: stress_test.py --detector detect_css.py --tag v3 [--hosts P001,P002,P003,P004,P005]
                      [--depths 0.6,0.45,0.3] [--seeds 1] [--z -2.5]
Inserts synthetic cSS into each host, runs the detector, and measures sensitivity, worst
rank, how many real non-lesion candidates (veins/artifacts) outrank the lesions, and AUC."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import argparse, os, re, shutil, subprocess, sys, numpy as np, pandas as pd, nibabel as nib
from scipy import ndimage as ndi
ap = argparse.ArgumentParser()
ap.add_argument("--detector", default="detect_css.py"); ap.add_argument("--tag", default="run")
ap.add_argument("--hosts", default="P001,P002,P003,P004,P005")
ap.add_argument("--depths", default="0.6,0.45,0.3"); ap.add_argument("--seeds", default="1")
ap.add_argument("--z", default="-2.5")
a = ap.parse_args()
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project")); S = f"{base}/scripts"
py = sys.executable
def run(*args):
    r = subprocess.run([py, *args], capture_output=True, text=True)
    if r.returncode: print(r.stdout[-500:], r.stderr[-1500:]); sys.exit("step failed: " + " ".join(args))
    return r.stdout
rows = []; feats = []
for h in a.hosts.split(","):
    clean = run(f"{S}/{a.detector}", h, a.z, "85", "3").split("\n")[0]
    nclean = int(re.search(r"-> (\d+) candidates", clean).group(1))
    shutil.copy(f"{base}/synthseg/{h}_seg.nii.gz", f"{base}/synthseg/{h}S_seg.nii.gz")
    for d in a.depths.split(","):
        for s in a.seeds.split(","):
            run(f"{S}/make_synthetic.py", h, "8", d, s)
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
                cdf["is_lesion"] = cdf.cand_id.isin(les).astype(int); cdf["host"] = h; cdf["depth"] = float(d)
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
              "surface_gradient", "surface_contact", "elongation", "score"]:
        if f in F.columns:
            L, N = F[F.is_lesion == 1][f], F[F.is_lesion == 0][f]
            print(f"  {f:18s} cSS median {L.median():7.2f} | others median {N.median():7.2f} | AUC {auc(L, N):.2f}")
```

## tests/phantom.py

```python
"""Synthetic phantom brain (ellipsoid, radial sulci, midline fissure, dark veins) for tests.
usage: phantom.py NAME SEED   (writes into $CSS_BASE)"""
import os, numpy as np, nibabel as nib, sys
B = os.environ["CSS_BASE"]
from scipy import ndimage as ndi
name, sd = sys.argv[1], int(sys.argv[2]); rng = np.random.default_rng(sd)
sh = (220, 220, 40); vox = (0.5, 0.5, 2.0)
aff = np.diag(list(vox) + [1.0])
x, y, zz = np.meshgrid(*[(np.arange(n) - n/2) * v for n, v in zip(sh, vox)], indexing="ij")
r = np.sqrt((x/50)**2 + (y/52)**2 + (zz/38)**2)          # ellipsoid "brain"
brain = r < 1
seg = np.zeros(sh, np.int32)
depth_mm = (1 - r) * 50                                    # approx depth from surface
left = x < 0; front = y > 0
ctxcode = np.where(left, 1000, 2000) + np.where(front, 28, 22)
seg[brain] = np.where(left, 2, 41)[brain]
seg[brain & (depth_mm < 3)] = ctxcode[brain & (depth_mm < 3)]
ang = np.arctan2(y, x)
I = np.zeros(sh, np.float32)
for a0 in np.arange(-np.pi, np.pi, np.pi/9):              # radial sulci
    d_line = np.abs(np.sin(ang - a0)) * np.sqrt(x**2 + y**2)
    sul = brain & (d_line < 0.8) & (depth_mm < 12) & (np.cos(ang - a0) > 0)
    bank = brain & (d_line < 3.8) & (depth_mm < 15) & (np.cos(ang - a0) > 0) & ~sul
    seg[bank] = ctxcode[bank]; seg[sul] = 24
seg[np.abs(x) < 0.7] = np.where(brain[np.abs(x) < 0.7], 24, 0)   # midline fissure
I[np.isin(seg, [2, 41])] = 230; I[seg >= 1000] = 200; I[seg == 24] = 255
I[brain] += rng.normal(0, 12, brain.sum())
# veins: dark tubes along some sulci, with a branch
for a0 in np.arange(-np.pi, np.pi, np.pi/4):
    d_line = np.abs(np.sin(ang - a0 + 0.17)) * np.sqrt(x**2 + y**2)
    v = brain & (d_line < 0.6) & (depth_mm > 2) & (depth_mm < 9) & (np.cos(ang - a0) > 0) & (np.abs(zz) < 14)
    I[v] = 70
I[~brain] = 0; I[brain] = np.maximum(I[brain], 1)
seg[~brain] = 0
nib.save(nib.Nifti1Image(I, aff), f"{B}/data/{name}_swi.nii")
nib.save(nib.Nifti1Image(seg, aff), f"{B}/synthseg/{name}_seg.nii.gz")
```

## tests/run_tests.sh

```bash
#!/bin/bash
# Self-test with a synthetic phantom brain - no patient data, no FreeSurfer needed.
# usage: bash tests/run_tests.sh      (from the package root, inside the 'css' conda env)
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
T=$(mktemp -d); export CSS_BASE=$T FREESURFER_HOME=$T/fs CSS_REVIEW_TEST=1
mkdir -p $T/data $T/synthseg $T/work $T/review $T/results $T/fs $T/scripts
cp $ROOT/scripts/*.py $T/scripts/
printf '0 Unknown 0 0 0 0\n2 Left-Cerebral-White-Matter 245 245 245 0\n24 CSF 60 60 60 0\n41 Right-Cerebral-White-Matter 245 245 245 0\n1022 ctx-lh-postcentral 220 20 20 0\n1028 ctx-lh-superiorfrontal 20 220 160 0\n2022 ctx-rh-postcentral 220 20 20 0\n2028 ctx-rh-superiorfrontal 20 220 160 0\n' > $T/fs/FreeSurferColorLUT.txt
S=$T/scripts
echo "== phantoms";            python $ROOT/tests/phantom.py PH1 1; python $ROOT/tests/phantom.py PH2 2
echo "== align";               python $S/align_seg.py PH1; python $S/align_seg.py PH2
echo "== detect (clean)";      python $S/detect_css.py PH1 | head -1
echo "== export";              python $S/export_review.py PH1 3
echo "== synthetic lesions";   python $S/make_synthetic.py PH1 8 0.6 1 | head -1
cp $T/synthseg/PH1_seg.nii.gz $T/synthseg/PH1S_seg.nii.gz; python $S/align_seg.py PH1S > /dev/null
python $S/detect_css.py PH1S | head -1
R=$(python $S/eval_synthetic.py PH1S | grep RESULT); echo "$R"
FOUND=$(echo "$R" | sed -E 's/.*found=([0-9]+).*/\1/')
[ "$FOUND" -ge 4 ] || { echo "FAIL: sensitivity too low ($FOUND/8)"; exit 1; }
echo "== scoring with truth";  python $S/score_css.py PH1S --truth > $T/score_test.log; grep -A2 "cSS multifocality" $T/score_test.log
echo "== interactive reviewer (headless test mode)"; python $S/review_css.py PH1S --top 8 --redo > $T/review_test.log; grep "cSS multifocality" $T/review_test.log
echo "== feature report";      python $S/feature_report.py | head -3
echo "== stress test (small)"; python $S/stress_test.py --hosts PH1,PH2 --depths 0.6 --seeds 1 --tag test | grep -A2 SUMMARY
echo; echo "ALL TESTS PASSED  (temp dir $T)"
```

