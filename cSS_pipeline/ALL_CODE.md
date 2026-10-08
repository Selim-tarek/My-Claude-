# All code — cSS pipeline v4

Generated from the files below; the files are authoritative.

## scripts/run_css.sh

```bash
#!/bin/zsh
# One patient, start to finish (up to review).
# usage: run_css.sh SUBJECT_ID [path/to/swi.nii(.gz)] [--nostrip] [--t1 T1.nii] [--flair FLAIR.nii] [--recon]
#   - with a SWI path: imports the scan; skull-strips it with SynthStrip unless --nostrip
#     (use --nostrip only for scans that are already skull-stripped)
#   - SWI = the processed SWI (or magnitude) series, NOT the minIP and NOT the phase map
#   - --t1 / --flair / --recon: v4 anatomy from T1 (and FLAIR), see prep_anat.sh
#     (T1 labels are re-used on later runs; SynthSeg on the SWI is then skipped)
#   - set NOVIEW=1 to skip opening freeview (batch use)
set -e
S=$1; shift 2>/dev/null || true
SRC=""; STRIP=1; ANAT=()
while [ $# -gt 0 ]; do
  case "$1" in
    --nostrip) STRIP=0; shift;;
    --t1|--flair) ANAT+=("$1" "$2"); shift 2;;
    --recon) ANAT+=("$1"); shift;;
    -*) echo "run_css.sh: unknown option $1"; exit 1;;
    *) SRC=$1; shift;;
  esac
done
B=${CSS_BASE:-$HOME/css_project}; SC=$B/scripts
if [ -z "$S" ]; then
  echo "usage: run_css.sh SUBJECT_ID [swi.nii] [--nostrip] [--t1 T1.nii] [--flair FLAIR.nii] [--recon]"; exit 1
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
if [ ${#ANAT[@]} -gt 0 ]; then
  echo "[0b] T1/FLAIR anatomy -> SWI grid"
  zsh $SC/prep_anat.sh $S "${ANAT[@]}"
fi
T1SEG=$B/work/${S}_t1seg_swispace.nii.gz
echo "[1/4] segmentation"
if [ -f $T1SEG ] && [ ! $B/data/${S}_swi.nii -nt $T1SEG ]; then
  echo "      using T1-derived labels ($T1SEG)"
# re-run when missing OR when the SWI was re-imported after the last segmentation (v4 fix:
# v3 silently reused a stale segmentation for a different scan with the same ID)
elif [ ! -f $B/synthseg/${S}_seg.nii.gz ] || [ $B/data/${S}_swi.nii -nt $B/synthseg/${S}_seg.nii.gz ]; then
  echo "      SynthSeg on the SWI (no T1 labels - consider --t1 for better anatomy)"
  mri_synthseg --i $B/data/${S}_swi.nii --o $B/synthseg/${S}_seg.nii.gz \
    --parc --robust --threads 8 > /dev/null
else
  echo "      SynthSeg on the SWI (already done, skipping)"
fi
echo "[2/4] aligning segmentation to SWI"
python $SC/align_seg.py $S
echo "[3/4] detecting cSS candidates"
python $SC/detect_css.py $S -2.5 85 3
if [ -z "$NOVIEW" ]; then
  echo "[4/4] opening viewer"
  freeview -v $B/data/${S}_swi.nii \
    $B/review/${S}_candidates.nii.gz:colormap=lut:opacity=0.6 > /dev/null 2>&1 &
fi
echo ""
echo "NEXT: python $SC/review_css.py $S          (interactive review of ALL candidates + score)"
echo "      or: python $SC/export_review.py $S 40 | pbcopy    (paste into the workbook)"
```

## scripts/prep_anat.sh

```bash
#!/bin/zsh
# Bring T1 (and FLAIR) anatomy onto the SWI grid (v4).
# usage: prep_anat.sh SUBJECT_ID --t1 T1.nii(.gz) [--flair FLAIR.nii(.gz)] [--recon]
#   --recon  full FreeSurfer recon-all on the T1 (about 1-3 h on an M3): gives Destrieux sulcal
#            labels -> true SULCAL multifocality scoring. Re-used if subjects/ID already finished.
#   default  SynthSeg --parc --robust on the T1 (minutes): sharper cortex/CSF than SynthSeg on SWI,
#            but no sulcal labels (scoring stays Euclidean).
# Needs data/ID_swi.nii (run_css.sh imports it). Writes, all on the SWI grid, in work/:
#   ID_t1seg_swispace.nii.gz    aseg + DK cortex labels  -> align_seg.py / detect_css.py
#   ID_a2009s_swispace.nii.gz   Destrieux labels         -> score_css.py (recon-all only)
#   ID_flair_swispace.nii.gz    FLAIR                    -> detect_css.py (acute cSAH check)
#   ID_swi2t1.lta, ID_flair2swi.lta   registrations: ALWAYS check them with the printed freeview line
# Use the 3-D T1 (MPRAGE / SPGR / BRAVO) and the 3-D or 2-D FLAIR. Do NOT use the minIP series.
set -e
S=$1; shift 2>/dev/null || true
T1=""; FL=""; RECON=0
while [ $# -gt 0 ]; do
  case "$1" in
    --t1) T1=$2; shift 2;;
    --flair) FL=$2; shift 2;;
    --recon) RECON=1; shift;;
    *) echo "prep_anat.sh: unknown option $1"; exit 1;;
  esac
done
if [ -z "$S" ] || [ -z "$T1" ]; then
  echo "usage: prep_anat.sh SUBJECT_ID --t1 T1.nii [--flair FLAIR.nii] [--recon]"; exit 1
fi
B=${CSS_BASE:-$HOME/css_project}; W=$B/work
# always the project folder (FreeSurfer's setup script points SUBJECTS_DIR at its own install dir);
# set CSS_SUBJECTS_DIR to use a different folder
export SUBJECTS_DIR=${CSS_SUBJECTS_DIR:-$B/subjects}
SWI=$B/data/${S}_swi.nii
[ -f $SWI ] || { echo "ERROR: $SWI not found - run run_css.sh $S <swi> first"; exit 1; }
mkdir -p $W $SUBJECTS_DIR

echo "[a] T1: import + skull-strip (SynthStrip)"
mri_convert "$T1" $B/data/${S}_t1.nii.gz > /dev/null
mri_synthstrip -i $B/data/${S}_t1.nii.gz -o $W/${S}_t1_brain.nii.gz > /dev/null

A2=""
if [ $RECON -eq 1 ] || [ -f $SUBJECTS_DIR/$S/mri/aparc.a2009s+aseg.mgz ]; then
  if [ ! -f $SUBJECTS_DIR/$S/mri/aparc.a2009s+aseg.mgz ]; then
    echo "[b] recon-all -all (1-3 h; log: $W/${S}_recon.log)"
    recon-all -s $S -i $B/data/${S}_t1.nii.gz -all -threads 8 > $W/${S}_recon.log 2>&1
  else
    echo "[b] recon-all already finished for $S - re-using it"
  fi
  LAB=$SUBJECTS_DIR/$S/mri/aparc+aseg.mgz
  A2=$SUBJECTS_DIR/$S/mri/aparc.a2009s+aseg.mgz
  echo "[c] register SWI -> T1 (bbregister, boundary-based, T2*-like contrast)"
  bbregister --s $S --mov $SWI --reg $W/${S}_swi2t1.lta --t2 --init-coreg > $W/${S}_bbreg.log 2>&1
  grep -i "min cost" $W/${S}_bbreg.log | tail -1 || true
else
  echo "[b] SynthSeg --parc --robust on the T1"
  mri_synthseg --i $B/data/${S}_t1.nii.gz --o $W/${S}_t1synthseg.nii.gz --parc --robust --threads 8 > /dev/null
  LAB=$W/${S}_t1synthseg.nii.gz
  echo "[c] register SWI -> T1 (mri_coreg, rigid, mutual information)"
  mri_coreg --mov $SWI --ref $W/${S}_t1_brain.nii.gz --reg $W/${S}_swi2t1.lta --dof 6 > $W/${S}_coreg.log 2>&1
fi

echo "[d] labels -> SWI grid (nearest neighbour)"
mri_vol2vol --mov $SWI --targ $LAB --lta $W/${S}_swi2t1.lta --inv --interp nearest \
  --o $W/${S}_t1seg_swispace.nii.gz > /dev/null
if [ -n "$A2" ]; then
  mri_vol2vol --mov $SWI --targ $A2 --lta $W/${S}_swi2t1.lta --inv --interp nearest \
    --o $W/${S}_a2009s_swispace.nii.gz > /dev/null
fi

if [ -n "$FL" ]; then
  echo "[e] FLAIR: import, skull-strip, register -> SWI, resample"
  mri_convert "$FL" $B/data/${S}_flair.nii.gz > /dev/null
  mri_synthstrip -i $B/data/${S}_flair.nii.gz -o $W/${S}_flair_brain.nii.gz > /dev/null
  mri_coreg --mov $W/${S}_flair_brain.nii.gz --ref $SWI --reg $W/${S}_flair2swi.lta --dof 6 > $W/${S}_flaircoreg.log 2>&1
  mri_vol2vol --mov $W/${S}_flair_brain.nii.gz --targ $SWI --lta $W/${S}_flair2swi.lta \
    --o $W/${S}_flair_swispace.nii.gz > /dev/null
fi

echo ""
echo "QC - labels must hug the SWI cortex, FLAIR must line up (toggle layers):"
echo "  freeview -v $SWI $W/${S}_t1seg_swispace.nii.gz:colormap=lut:opacity=0.3" \
     "$( [ -n "$FL" ] && echo "$W/${S}_flair_swispace.nii.gz:visible=0" )" \
     "$( [ -n "$A2" ] && echo "$W/${S}_a2009s_swispace.nii.gz:colormap=lut:opacity=0.3:visible=0" )"
echo "If misaligned: delete work/${S}_t1seg_swispace.nii.gz and the pipeline falls back to SWI SynthSeg."
```

## scripts/run_all.sh

```bash
#!/bin/zsh
# Batch: run several patients without opening the viewer.
# usage: run_all.sh P006 P007 P008        (expects data/<ID>_swi.nii already imported)
# Full output of each subject goes to $CSS_BASE/logs/<ID>_run.log; the detector summary is printed.
# (v4 fix: v3 piped into `head -1`, which closed the pipe and killed run_css.sh after SynthSeg,
#  so detection never ran in batch mode.)
B=${CSS_BASE:-$HOME/css_project}
mkdir -p $B/logs
for s in "$@"; do
  echo "======== $s"
  if NOVIEW=1 zsh $B/scripts/run_css.sh $s > $B/logs/${s}_run.log 2>&1; then
    grep -m1 "raw components" $B/logs/${s}_run.log
  else
    echo "FAILED - see $B/logs/${s}_run.log"; tail -5 $B/logs/${s}_run.log
  fi
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

## scripts/css_common.py

```python
"""Shared helpers for the cSS pipeline scripts (imported from the same folder)."""
import os, numpy as np, pandas as pd

FP = ["vox_i", "vox_j", "vox_k", "volume_mm3"]      # candidate fingerprint stored with each call


def base():
    return os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))


def load_lut():
    lut = {}
    p = os.path.join(os.environ.get("FREESURFER_HOME", ""), "FreeSurferColorLUT.txt")
    if os.path.exists(p):
        for line in open(p):
            t = line.split()
            if len(t) > 2 and t[0].isdigit():
                lut[int(t[0])] = t[1]
    return lut


def match_calls(cand, calls, tol_vox=2.0, tol_vol=0.15):
    """Re-attach reader calls to the CURRENT candidate list by fingerprint (centroid within
    tol_vox voxels, volume within tol_vol), so re-running the detector can never move a call
    onto a different lesion. Returns ({cand_id: call}, n_dropped), or (None, n) for a legacy
    calls file that has no fingerprint columns."""
    if not all(c in calls.columns for c in FP):
        return None, len(calls)
    out, dropped = {}, 0
    C = cand[["cand_id"] + FP].to_numpy(float)
    for _, r in calls.iterrows():
        d = np.sqrt(((C[:, 1:4] - r[FP[:3]].to_numpy(float)) ** 2).sum(1))
        dv = np.abs(C[:, 4] - float(r.volume_mm3)) / max(float(r.volume_mm3), 1.0)
        ok = np.where((d <= tol_vox) & (dv <= tol_vol))[0]
        if len(ok):
            out[int(C[ok[d[ok].argmin()], 0])] = r.call
        else:
            dropped += 1
    return out, dropped


def load_calls(b, subj, cand):
    """Reader calls for subj mapped onto the current candidates. Returns (dict, warning)."""
    p = f"{b}/review/{subj}_calls.csv"
    if not os.path.exists(p) or not len(cand):
        return {}, ""
    calls = pd.read_csv(p)
    if not len(calls):
        return {}, ""
    m, dropped = match_calls(cand, calls)
    if m is None:
        # v3 file (cand_id only): trust it only if the detector has not been re-run since.
        # candidates.nii.gz is written by detect_css.py only (mark_css.py rewrites the csv).
        if os.path.getmtime(f"{b}/review/{subj}_candidates.nii.gz") > os.path.getmtime(p):
            return {}, (f"{subj}: old-format calls IGNORED - candidates were re-detected after "
                        f"they were made, so cand_ids may point to different lesions")
        valid = set(cand.cand_id.astype(int))
        return {int(c): v for c, v in zip(calls.cand_id, calls.call) if int(c) in valid}, ""
    return m, (f"{subj}: {dropped} earlier call(s) match no current candidate (detector re-run) "
               f"- dropped" if dropped else "")


def save_calls(b, subj, cand, calls):
    c = cand.set_index("cand_id")
    rows = [dict(cand_id=k, call=v, **{f: c.loc[k, f] for f in FP}) for k, v in sorted(calls.items())]
    pd.DataFrame(rows, columns=["cand_id", "call"] + FP).to_csv(f"{b}/review/{subj}_calls.csv", index=False)
```

## scripts/align_seg.py

```python
"""Labels on the SWI grid -> work/ID_seg_swispace.nii.gz (masked to SWI coverage).
Source, best first:
  work/ID_t1seg_swispace.nii.gz  T1-derived labels already in SWI space (prep_anat.sh, v4)
  synthseg/ID_seg.nii.gz         SynthSeg run on the SWI itself (v3 behaviour)"""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib
from nibabel.processing import resample_from_to
from scipy import ndimage as ndi
subj = sys.argv[1]
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
swi = nib.load(f"{base}/data/{subj}_swi.nii")
t1p = f"{base}/work/{subj}_t1seg_swispace.nii.gz"
src = t1p if os.path.exists(t1p) and os.path.getmtime(t1p) >= os.path.getmtime(f"{base}/data/{subj}_swi.nii") \
    else f"{base}/synthseg/{subj}_seg.nii.gz"
seg = nib.load(src)
seg_swi = resample_from_to(seg, swi, order=0).get_fdata().astype(np.int32)
n_before = (seg_swi > 0).sum()
seg_swi[~ndi.binary_fill_holes(swi.get_fdata() > 0)] = 0
out = nib.Nifti1Image(seg_swi, swi.affine, swi.header); out.set_data_dtype(np.int32)
nib.save(out, f"{base}/work/{subj}_seg_swispace.nii.gz")
print(f"{subj}: labels from {os.path.basename(src)}; labelled voxels before mask {n_before:,}, "
      f"after mask {(seg_swi>0).sum():,}")
```

## scripts/detect_css.py

```python
#!/usr/bin/env python3
"""cSS candidate detector v4 - definition-driven.
usage: detect_css.py SUBJECT [z_thr=-2.5] [ridge_pct=85] [rim_mm=3]
env:   CSS_RANK=v3 ranks with the old v3 score (both scores are always written)

cSS (STRIVE-2 / Charidimou): curvilinear hypointensity in the SUBPIAL cortex / subarachnoid
space that follows the gyral surface, often on both banks of a sulcus ("tram-track"),
supratentorial, remote from lobar ICH, no FLAIR hyperintensity (that would be acute cSAH).
Each criterion is measured per candidate:

  definition                      column(s)                      cSS            vein / mimic
  outlines the cortical surface   pial_dist_mm, bank_frac        ~0 mm, high    sulcal centre (>0)
  curvilinear sheet, not a tube   tube_ratio (3-D Hessian)       low (plate)    high (tube)
  follows the gyral contour       surface_alignment              high (~1)      lower
  tram-track (both banks)         tram_frac                      >0             ~0
  not part of the venous tree     vein_tree_mm (recorded only)   short          long network
  asymmetric                      mirror_dark_frac               low            high (normal veins)
  remote from ICH                 ich_dist_mm, near_ich_suggest  far            <=5 mm
  supratentorial                  infratentorial                 0              1 (classical SS)
  chronic, not acute cSAH         flair_csf_z, flair_bright      FLAIR not bright (hint only:
                                                                 chronic cSS can show mild FLAIR signal)
  mimic: cortical vein thrombosis flair_ctx_z, flair_ctx_bright  FLAIR-bright cortex nearby -> check
                                  (only if work/ID_flair_swispace.nii.gz exists)

Candidate generation (unchanged idea from v3, tightened):
  - search zone = pial band: cortex + subarachnoid CSF, -3.5..+4 mm from the pial surface, minus
    white matter / deep nuclei (v3 dilated cortex into WM, which is darker than GM on SWI)
  - darkness z vs a LOCAL cortex reference that now EXCLUDES locally dark voxels (2-pass), so
    extensive cSS no longer darkens its own reference
  - 2-D Sato ridge + hysteresis, per hemisphere; >=25 mm3, >=2 slices, elongation >=2.0
Labels: work/ID_seg_swispace.nii.gz (from T1 via prep_anat.sh when available, else SWI SynthSeg).
Lobar ICH mask: work/ID_ich.nii.gz if you drew one (freeview), else found automatically.
"""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
from skimage.filters import sato
from skimage.measure import regionprops
from skimage.morphology import skeletonize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, load_lut

subj   = sys.argv[1]
z_thr  = float(sys.argv[2]) if len(sys.argv) > 2 else -2.5
r_pct  = float(sys.argv[3]) if len(sys.argv) > 3 else 85
RIM_MM = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
RANK   = os.environ.get("CSS_RANK", "v4")
Z_LOW, R_PCT_LOW = -2.0, 80           # hysteresis (weak) thresholds
GROW_Z = -1.5                          # 'dark' map used for region growing at scoring
MIN_MM3, MIN_SLICES, MIN_ELONG = 25.0, 2, 2.0   # v4: 2.0 (tram-track = two parallel plates, ~2.4;
                                                # round microbleeds stay ~1-1.5)
PIAL_IN_MM, PIAL_OUT_MM = 3.5, 4.0     # pial band: cortex thickness + blooming / subarachnoid CSF
HESS_SIGMAS_MM = (0.8, 1.6)            # 3-D shape analysis scales
TRAM_MAX_MM = 10.0                     # look for the opposite bank up to this far across a sulcus
ICH_MIN_MM3, ICH_MIN_RADIUS_MM = 500.0, 2.5   # automatic lobar ICH: dark, compact (inscribed
                                               # radius >= 2.5 mm, i.e. not a vein/sheet) blob
ARTIFACT_ZONES = ("orbitofrontal", "entorhinal", "temporalpole", "parahippocampal",
                  "fusiform", "inferiortemporal")
LH = [2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 17, 18, 26, 28]
RH = [41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60]
WM = [2, 41, 77, 251, 252, 253, 254, 255]
DEEP = [4, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 26, 28, 31, 43, 44, 49, 50, 51, 52, 53, 54, 58, 60, 63]
INFRA = [7, 8, 46, 47, 16, 15]
CSF = [24]

base = _base()
swi = nib.load(f"{base}/data/{subj}_swi.nii")
I   = swi.get_fdata().astype(np.float32)
seg = nib.load(f"{base}/work/{subj}_seg_swispace.nii.gz").get_fdata().astype(np.int32)
vox = tuple(float(v) for v in swi.header.get_zooms()[:3])
s2d = np.ones((3, 3, 1), bool); s3 = np.ones((3, 3, 3), bool)
mm2vox = lambda mm: tuple(mm / v for v in vox)

brain  = ndi.binary_fill_holes(I > 0)
cortex = (seg >= 1000) | (seg == 3) | (seg == 42)
it = max(1, int(round(2.0 / vox[0])))
edge  = brain & ~ndi.binary_erosion(brain, structure=s2d, iterations=it)
surface = (~brain) | np.isin(seg, CSF)
near_surface = ndi.binary_dilation(surface, structure=s2d, iterations=it)

# ---- 0. signed distance to the pial surface (mm; <0 inside tissue, >0 in CSF) + normals
tissue = ndi.binary_fill_holes(brain & (seg > 0) & ~np.isin(seg, CSF))    # fills ventricles
pial_sd = (ndi.distance_transform_edt(~tissue, sampling=vox)
           - ndi.distance_transform_edt(tissue, sampling=vox)).astype(np.float32)
_g = np.gradient(ndi.gaussian_filter(pial_sd, mm2vox(1.0)), *vox)
_gn = np.sqrt(sum(g * g for g in _g)) + 1e-6
normal = [(g / _gn).astype(np.float32) for g in _g]          # points from tissue into CSF
del _g, _gn

# ---- 1. robust normalisation, local cortex reference (2-pass, dark voxels excluded)
hi = np.percentile(I[brain], 99.5)
Ic = np.minimum(I, hi)
c = Ic[cortex & (I > 0)]
mu = float(np.median(c)); sd = float(1.4826 * np.median(np.abs(c - mu)) + 1e-6)
_sig = mm2vox(10.0)
def local_ref(ref):
    w = ndi.gaussian_filter(ref.astype(np.float32), _sig)
    m = ndi.gaussian_filter((Ic * ref).astype(np.float32), _sig)
    return np.where(w > 0.05, m / np.maximum(w, 1e-6), mu)
z = (Ic - local_ref(cortex)) / sd
# pass 2: locally dark cortex (lesion candidates) no longer pulls its own reference down;
# regionally dark NORMAL cortex (iron-rich motor cortex) is not locally dark, so it stays in
z = (Ic - local_ref(cortex & (z > Z_LOW))) / sd
z_glob = (Ic - mu) / sd

# ---- 2. search zone: pial band minus skull-strip rim, white matter and deep nuclei
dist_mm = ndi.distance_transform_edt(brain, sampling=vox)
rim = brain & (dist_mm < RIM_MM)
If = Ic.copy(); If[~brain] = mu
shell = brain & (pial_sd >= -PIAL_IN_MM) & (pial_sd <= PIAL_OUT_MM) & ~np.isin(seg, WM + DEEP)
shell |= edge & ~np.isin(seg, WM + DEEP)
zone = shell & ~rim

# ---- 3. curvilinear dark structures, slice by slice (candidate generation)
R = np.zeros_like(I)
for k in range(I.shape[2]):
    if I[:, :, k].max() > 0:
        R[:, :, k] = sato(If[:, :, k], sigmas=[1, 2, 3], black_ridges=True)
r_hi = np.percentile(R[shell], r_pct); r_lo = np.percentile(R[shell], R_PCT_LOW)
strong = zone & (z < z_thr) & (R > r_hi)
weak   = zone & (z < Z_LOW) & (R > r_lo)          # strong is a subset of weak

# ---- 4. hemisphere map; hysteresis labelling separately in each hemisphere
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
del R, wl, weak, strong

os.makedirs(f"{base}/work", exist_ok=True)
nib.save(nib.Nifti1Image((zone & (z < GROW_Z)).astype(np.uint8), swi.affine),
         f"{base}/work/{subj}_dark.nii.gz")

# ---- 5. 3-D Hessian shape (multi-scale) at dark voxels: plate (cSS) vs tube (vein)
need = brain & (z < Z_LOW) & (pial_sd > -PIAL_IN_MM)      # candidates + venous network voxels
need |= lab > 0
idx = np.nonzero(need)
best = np.full(len(idx[0]), -np.inf, np.float32)
eig = np.zeros((len(idx[0]), 3), np.float32); evec3 = np.zeros((len(idx[0]), 3), np.float32)
for s_mm in HESS_SIGMAS_MM:
    sg = mm2vox(s_mm); H = np.zeros((len(idx[0]), 3, 3), np.float32)
    for a in range(3):
        for b in range(a, 3):
            order = [0, 0, 0]; order[a] += 1; order[b] += 1
            h = ndi.gaussian_filter(Ic, sg, order=order)[idx] / (vox[a] * vox[b]) * s_mm ** 2
            H[:, a, b] = h; H[:, b, a] = h
    w, v = np.linalg.eigh(H)                              # ascending signed eigenvalues
    o = np.argsort(np.abs(w), axis=1)                     # sort by magnitude |l1|<=|l2|<=|l3|
    w = np.take_along_axis(w, o, 1); v3 = v[np.arange(len(v)), :, o[:, 2]]   # eigvec of l3
    upd = w[:, 2] > best                                  # dark structure: largest curvature > 0
    best[upd] = w[upd, 2]; eig[upd] = w[upd]; evec3[upd] = v3[upd]
    del H, w, v
l2, l3 = np.abs(eig[:, 1]), eig[:, 2]
tube_v = np.where(l3 > 0, l2 / np.maximum(np.abs(l3), 1e-6), np.nan).astype(np.float32)
TUBE = np.full(I.shape, np.nan, np.float32); TUBE[idx] = tube_v
ALIGN = np.full(I.shape, np.nan, np.float32)
ALIGN[idx] = np.abs(sum(evec3[:, a] * normal[a][idx] for a in range(3)))
del eig, evec3, best, tube_v

# ---- 6. venous network: dark tubular voxels in the subarachnoid space / cortex
vnet = need & (z < Z_LOW) & (TUBE > 0.5) & (pial_sd > -1.5)
vl, nv = ndi.label(vnet | (lab > 0), structure=s3)
vext = np.zeros(nv + 1, np.float32)
for n, sl in enumerate(ndi.find_objects(vl), 1):
    if sl is not None:
        vext[n] = np.sqrt(sum(((s.stop - s.start) * v) ** 2 for s, v in zip(sl, vox)))

# ---- 7. mid-sagittal plane (hemisphere boundary) for the mirror check
bd = brain & hemiL & ndi.binary_dilation(~hemiL, structure=s3)
P = np.argwhere(bd).astype(np.float32) * vox
if len(P) > 100:
    pc = P.mean(0); mvec = np.linalg.svd(P[:: max(1, len(P) // 20000)] - pc, full_matrices=False)[2][2]
else:
    pc = np.array(I.shape) * vox / 2; mvec = np.array([1.0, 0, 0])
zmin3 = ndi.minimum_filter(z, size=3)
def mirror_dark(pts):
    p = pts * vox; q = p - 2 * ((p - pc) @ mvec)[:, None] * mvec[None, :]
    qi = np.clip(np.round(q / vox).astype(int), 0, np.array(I.shape) - 1)
    return float((zmin3[qi[:, 0], qi[:, 1], qi[:, 2]] < Z_LOW).mean())

# ---- 8. lobar ICH: reader-drawn mask, else dark compact blobs in lobar tissue
ich_p = f"{base}/work/{subj}_ich.nii.gz"
lut = load_lut()
ART_CODES = [c for c, n in lut.items() if any(a in n for a in ARTIFACT_ZONES)]
if os.path.exists(ich_p):
    ich = nib.load(ich_p).get_fdata() > 0; ich_src = "drawn"
else:
    wm_mu = float(np.median(Ic[np.isin(seg, [2, 41])])) if np.isin(seg, [2, 41]).any() else mu
    lobar = np.isin(seg, [2, 41]) | cortex
    cand_ich = lobar & ((Ic - wm_mu) / sd < -3.0) & (pial_sd < -1.0)
    il, ni = ndi.label(cand_ich, structure=s3)
    ich = np.zeros_like(brain)
    for n, sl in enumerate(ndi.find_objects(il), 1):
        if sl is None: continue
        m = il[sl] == n
        if m.sum() * np.prod(vox) < ICH_MIN_MM3: continue
        if ndi.distance_transform_edt(np.pad(m, 1), sampling=vox).max() < ICH_MIN_RADIUS_MM: continue
        # a lobar haematoma lies in the tissue under the cortex. On P006 the auto-finder caught
        # skull-base susceptibility artifact (orbitofrontal) and the sagittal sinus / vertex veins:
        # reject blobs that are superficial, midline, at the skull-strip edge or in artifact zones
        if np.median(pial_sd[sl][m]) > -3.0: continue
        if np.median(dmid[sl][m]) < 5.0: continue
        if rim[sl][m].mean() > 0.2: continue
        if np.isin(seg[sl][m], ART_CODES).mean() > 0.3: continue
        ich[sl] |= m
    ich_src = "auto"
ich_dist = ndi.distance_transform_edt(~ich, sampling=vox) if ich.any() else np.full(I.shape, np.inf)

# ---- 9. optional FLAIR (acute convexity SAH is FLAIR-bright in the sulci; chronic cSS is not)
fl_p = f"{base}/work/{subj}_flair_swispace.nii.gz"
FL = None
if os.path.exists(fl_p):
    FL = nib.load(fl_p).get_fdata().astype(np.float32)
    csf_ref = FL[np.isin(seg, CSF) & brain & (FL > 0)]
    tis_ref = FL[tissue & brain & (FL > 0)]
    ctx_ref = FL[cortex & brain & (FL > 0)]
    if csf_ref.size and tis_ref.size and ctx_ref.size:
        f_mu = float(np.median(csf_ref)); f_sd = float(1.4826 * np.median(np.abs(tis_ref - np.median(tis_ref))) + 1e-6)
        fc_mu = float(np.median(ctx_ref))
    else:
        FL = None

# ---- 10. features per candidate
cortex_dist = ndi.distance_transform_edt(~cortex, sampling=vox)      # 0 inside cortex
pial_depth = ndi.distance_transform_edt(cortex | np.isin(seg, [2, 41]), sampling=vox)
step = 0.4; T_ = np.arange(step, TRAM_MAX_MM + step, step, dtype=np.float32)

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

def tram_track(pts):
    """From each bank voxel, walk along the surface normal across the sulcal CSF to the opposite
    bank. tram = opposite bank also dark AND the sulcal centre between them brighter than both
    (two separate dark lines = tram-track; a vein in the middle gives one central line)."""
    sdv = pial_sd[tuple(pts.T)]
    pts = pts[(sdv > -2.5) & (sdv < 1.0)]
    if len(pts) < 5: return np.nan, 0
    if len(pts) > 400: pts = pts[np.linspace(0, len(pts) - 1, 400).astype(int)]
    n = np.stack([normal[a][tuple(pts.T)] for a in range(3)], 1)
    q = pts[:, None, :] + T_[None, :, None] * n[:, None, :] / np.array(vox)   # voxel coords
    co = q.reshape(-1, 3).T
    sdr = ndi.map_coordinates(pial_sd, co, order=1, mode="nearest").reshape(q.shape[:2])
    zr = ndi.map_coordinates(z, co, order=1, mode="nearest").reshape(q.shape[:2])
    zs = z[tuple(pts.T)]
    hits = tram = 0
    for i in range(len(pts)):
        out = np.nonzero(sdr[i] > 0.3)[0]
        if not len(out): continue
        back = np.nonzero((sdr[i] <= 0) & (np.arange(len(T_)) > out[0]))[0]
        if not len(back): continue
        j = back[0]; hits += 1
        zo = zr[i, max(j - 2, 0): j + 5].min()
        zc_ = np.median(zr[i, out[0]: j][np.argsort(-sdr[i, out[0]: j])[:3]])   # sulcal midline
        if zo < Z_LOW and zc_ > max(zs[i], zo) + 1.0: tram += 1
    return (tram / hits if hits >= 5 else np.nan), hits

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
    seg_nb = seg[pad][ndi.binary_dilation(mm, iterations=2)]
    nb = seg_nb[seg_nb >= 1000]
    # infratentorial only when clearly cerebellar/brainstem with little cerebral cortex nearby, so
    # occipital / inferior temporal cSS next to the tentorium is not dropped from the 0-4 score
    infra = bool(np.isin(seg_nb, INFRA).mean() > 0.5 and nb.size < 0.2 * seg_nb.size)
    if nb.size == 0 and not infra: continue
    # hemisphere = majority of the candidate's own voxels in the hemisphere map (v4 fix: v3 counted
    # own L/R labels, but cSS lies in sulcal CSF where both counts are 0, so medial candidates
    # could take the region - and the score - of the opposite hemisphere)
    is_left = bool(hemiL[pad][mm].mean() >= 0.5)
    if nb.size:
        side = nb[nb < 2000] if is_left else nb[nb >= 2000]
        if side.size: nb = side
        code = int(np.bincount(nb).argmax())
    else:
        code = int(np.bincount(seg_nb[np.isin(seg_nb, INFRA)]).argmax())
    region = lut.get(code, str(code))
    pts = np.argwhere(mm) + [s.start for s in pad]
    contact   = float(near_surface[pad][mm].mean())
    depth     = float(-z[pad][mm].mean())
    cfrac     = float(cortex[pad][mm].mean())
    cdist     = float(np.median(cortex_dist[pad][mm]))
    bdens     = branches_per_10mm(mm)
    nbh = ndi.binary_dilation(mm, structure=s2d, iterations=4) & cortex[pad]
    dp = pial_depth[pad]; zz = z[pad]
    outer = nbh & (dp <= 1.0); inner = nbh & (dp > 1.0) & (dp <= 2.5)
    sgrad = float(zz[inner].mean() - zz[outer].mean()) if outer.sum() >= 5 and inner.sum() >= 5 else float('nan')
    # v4 definition features
    psd = pial_sd[pad][mm]
    pdist = float(np.median(psd))
    bankf = float(((psd > -2.5) & (psd < 1.0)).mean())
    tub = float(np.nanmedian(TUBE[pad][mm])) if np.isfinite(TUBE[pad][mm]).any() else np.nan
    aln = float(np.nanmedian(ALIGN[pad][mm])) if np.isfinite(ALIGN[pad][mm]).any() else np.nan
    tram, _ = tram_track(pts)
    comps = np.unique(vl[pad][ndi.binary_dilation(mm, structure=s3)]); comps = comps[comps > 0]
    vtree = float(vext[comps].max()) if len(comps) else 0.0
    mdark = mirror_dark(pts)
    idist = float(ich_dist[pad][mm].min())
    fz = fcz = np.nan
    if FL is not None:
        near3 = ndi.binary_dilation(mm, structure=s3, iterations=max(1, int(round(3 / vox[0])))) & (FL[pad] > 0)
        ring = near3 & (pial_sd[pad] > 0.5)                       # sulcal CSF: acute convexity SAH
        if ring.sum() >= 5: fz = float(np.median((FL[pad][ring] - f_mu) / f_sd))
        ctxn = near3 & cortex[pad]                                # cortex: oedema of cortical vein thrombosis
        if ctxn.sum() >= 5: fcz = float(np.median((FL[pad][ctxn] - fc_mu) / f_sd))
    artz      = any(a in region for a in ARTIFACT_ZONES)
    midz      = bool(np.median(dmid[pad][mm]) < 5.0)
    # off-surface tube. vein_tree_mm is recorded only: on real SWI (P006) the dark tubular network
    # along the cortex merges into one brain-wide component (255 mm for most candidates)
    vein_like = bool(pdist > 0.3 and (np.isnan(tub) or tub > 0.4))
    longs     = bool(n_sl * vox[2] > 40.0)    # > 40 mm top-to-bottom (sinus / large veins)
    # v3 score (kept for comparison)
    s3_ = depth * (0.25 + contact) * np.log1p(vol)
    for flag, f in ((artz, 0.5), (midz, 0.5), (cfrac < 0.35, 0.6), (longs, 0.5)):
        if flag: s3_ *= f
    # v4 score: darkness x how well each part of the definition is met (transparent, untrained)
    on_surf = float(np.exp(-((pdist + 0.5) / 1.5) ** 2))          # subpial: ~ -0.5 mm
    sheet = 1.0 - (tub if np.isfinite(tub) else 0.5)
    # sqrt: darkness is damped - on real data (P006) veins were DARKER than cSS (AUC 0.37)
    s4 = np.sqrt(max(depth, 0.0)) * (0.25 + on_surf) * (0.25 + sheet) * (0.5 + (aln if np.isfinite(aln) else 0.5)) \
         * (1.0 + (tram if np.isfinite(tram) else 0.0))
    for flag, f in ((artz, 0.5), (midz, 0.7), (mdark > 0.6, 0.8),
                    (longs, 0.5), (infra, 0.3)):
        if flag: s4 *= f
    c0 = pts.mean(0)
    rows.append(dict(hemi="L" if is_left else "R", region=region, label=code,
                     volume_mm3=round(vol, 1), n_slices=n_sl, elongation=round(elong, 1),
                     surface_contact=round(contact, 2), darkness_z=round(depth, 1),
                     artifact_zone=int(artz), midline_zone=int(midz),
                     score=round(s4 if RANK == "v4" else s3_, 2),
                     vox_i=int(round(c0[0])), vox_j=int(round(c0[1])), vox_k=int(round(c0[2])),
                     cortex_frac=round(cfrac, 2), cortex_dist_mm=round(cdist, 2),
                     branch_per10mm=round(bdens, 2), surface_gradient=round(sgrad, 2),
                     vein_like=int(vein_like), long_structure=int(longs),
                     pial_dist_mm=round(pdist, 2), bank_frac=round(bankf, 2),
                     tube_ratio=round(tub, 2), surface_alignment=round(aln, 2),
                     tram_frac=round(tram, 2), vein_tree_mm=round(vtree, 1),
                     mirror_dark_frac=round(mdark, 2), ich_dist_mm=round(min(idist, 999.0), 1),
                     near_ich_suggest=int(idist <= 5.0), infratentorial=int(infra),
                     flair_csf_z=round(fz, 2), flair_bright=int(np.isfinite(fz) and fz > 3.0),
                     flair_ctx_z=round(fcz, 2), flair_ctx_bright=int(np.isfinite(fcz) and fcz > 3.0),
                     score_v3=round(s3_, 2), score_v4=round(s4, 2), _lab=p.label))

# CSV columns up to "long_structure" are consumed by export_review.py / the workbook: keep them
cols = ["cand_id", "hemi", "region", "label", "volume_mm3", "n_slices", "elongation",
        "surface_contact", "darkness_z", "artifact_zone", "midline_zone", "score",
        "vox_i", "vox_j", "vox_k", "cortex_frac", "cortex_dist_mm", "branch_per10mm", "surface_gradient",
        "vein_like", "long_structure",
        "pial_dist_mm", "bank_frac", "tube_ratio", "surface_alignment", "tram_frac", "vein_tree_mm",
        "mirror_dark_frac", "ich_dist_mm", "near_ich_suggest", "infratentorial", "flair_csf_z",
        "flair_bright", "flair_ctx_z", "flair_ctx_bright", "score_v3", "score_v4", "accept"]
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
nib.save(nib.Nifti1Image(ich.astype(np.uint8), swi.affine), f"{base}/work/{subj}_ich_used.nii.gz")
cnt = lambda c: int(df[c].astype(int).sum()) if len(df) else 0
print(f"{subj}: {n_raw} raw components -> {len(df)} candidates "
      f"({cnt('artifact_zone')} artifact zone, {cnt('midline_zone')} midline, {cnt('vein_like')} vein-like, "
      f"{cnt('near_ich_suggest')} near ICH [{ich_src}], {cnt('infratentorial')} infratentorial)  rank={RANK}"
      + ("  FLAIR used" if FL is not None else "") + "\n")
if len(df):
    print(df.head(15)[["cand_id", "hemi", "region", "volume_mm3", "darkness_z", "pial_dist_mm",
                       "tube_ratio", "surface_alignment", "tram_frac", "vein_tree_mm", "score",
                       "vox_i", "vox_j", "vox_k"]].to_string(index=False))
```

## scripts/review_css.py

```python
#!/usr/bin/env python3
"""Interactive cSS review: shows each suspected candidate, you agree or disagree.
usage: review_css.py SUBJECT [--top N] [--redo]
  default: ALL candidates (v4). With --top N, candidates below N are recorded as NOT reviewed and
  the score is reported with that caveat (v3 silently scored them as "not cSS").

Keys (or click the buttons):
  y = cSS (agree)      v = vein      o = normal cortex      a = artifact
  i = near ICH         u = unsure    b = back one            [ / ] = move slice down/up
  q = finish and score
Decisions are saved after every key press, so you can quit and resume later.
Left column: whole slice (top) and an 8 mm minIP slab of the zoom window (bottom) computed from the
SWI itself - veins become continuous branching tubes there, cSS stays a band along the cortex.
(Do not feed the scanner minIP series into the pipeline; this panel is only for reading.)
Title: the v4 definition features - on-surface distance, plate/tube shape, tram-track, vein tree,
nearby ICH, FLAIR - as a reading aid; the call is always yours.
At the end, accepted candidates are scored automatically (score_css.py)."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, subprocess, numpy as np, nibabel as nib, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import load_calls, save_calls
import matplotlib
if os.environ.get("CSS_REVIEW_TEST"): matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.patches import Rectangle

subj = sys.argv[1]
top = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else None
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
if not redo:
    # v4 fix: calls are matched by candidate fingerprint (centroid + volume), not by rank, so a
    # re-run of the detector can no longer attach an old call to a different lesion
    calls, warn = load_calls(base, subj, df)
    if warn: print("WARNING:", warn)
KEYS = {"y": "cSS", "v": "Vein", "o": "Normal", "a": "Artifact", "i": "Near ICH", "u": "Unsure"}
COL = {"cSS": "#2e7d32", "Vein": "#1565c0", "Normal": "#6d6d6d", "Artifact": "#ef6c00",
       "Near ICH": "#8e24aa", "Unsure": "#c9a400"}
HW = int(round(30 / vox[0]))          # 60 mm zoom window
SLAB = max(1, int(round(4 / vox[2]))) # minIP slab half-thickness (~8 mm total)
def num(v):
    try: return float(v)
    except (TypeError, ValueError): return float("nan")
flag = lambda v: num(v) > 0          # NaN / missing column -> False

class Reviewer:
    def __init__(self):
        self.pos = next((n for n, c in enumerate(ids) if c not in calls), 0)
        self.shift = 0
        self.fig = plt.figure(figsize=(14, 8.2))
        gs = self.fig.add_gridspec(2, 4, left=0.02, right=0.98, top=0.85, bottom=0.14, wspace=0.05, hspace=0.12)
        self.ov = self.fig.add_subplot(gs[0, 0]); self.mip = self.fig.add_subplot(gs[1, 0])
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
        # minIP slab (8 mm) of the zoom window, suspect outlined on the centre slice
        k0, k1 = max(kc - SLAB, 0), min(kc + SLAB + 1, I.shape[2])
        slab = np.where(I[i0:i1, j0:j1, k0:k1] > 0, I[i0:i1, j0:j1, k0:k1], hi).min(axis=2)
        self.mip.clear(); self.mip.imshow(slab.T, cmap="gray", origin="lower", vmin=lo, vmax=hi)
        if m[i0:i1, j0:j1, k0:k1].any():
            self.mip.contour(m[i0:i1, j0:j1, k0:k1].any(axis=2).T.astype(float), levels=[0.5], colors="red", linewidths=0.8)
        self.mip.set_title(f"minIP {int(round((k1 - k0) * vox[2]))} mm slab (veins = tubes)", fontsize=9); self.mip.axis("off")
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
        flags = [n for n, f in (("artifact zone", "artifact_zone"), ("midline", "midline_zone"), ("vein-like", "vein_like"),
                                ("INFRATENTORIAL", "infratentorial"), ("FLAIR-bright CSF: acute cSAH?", "flair_bright"),
                                ("FLAIR-bright cortex: cortical vein thrombosis?", "flair_ctx_bright")) if flag(r.get(f, 0))]
        if flag(r.get("near_ich_suggest", 0)):
            flags.append(f"ICH {num(r.get('ich_dist_mm')):.0f} mm away - press i if contiguous")
        feat = ""
        if "pial_dist_mm" in r:
            shape = "plate" if num(r.tube_ratio) < 0.4 else ("tube" if num(r.tube_ratio) > 0.5 else "mixed")
            tram = num(r.tram_frac)
            feat = (f"\npial {num(r.pial_dist_mm):+.1f} mm   shape {shape} ({num(r.tube_ratio):.2f})   "
                    f"follows surface {num(r.surface_alignment):.2f}   "
                    f"tram-track {'n/a' if np.isnan(tram) else f'{tram:.0%}'}   vein tree {num(r.vein_tree_mm):.0f} mm   "
                    f"mirror dark {num(r.mirror_dark_frac):.0%}")
        done = sum(1 for c in ids if c in calls)
        prev = calls.get(cid)
        self.fig.suptitle(f"{subj}   candidate #{cid}  ({self.pos + 1}/{len(ids)}, {done} decided)   —   "
                          f"{r.hemi} {r.region}   {r.volume_mm3} mm³, {int(r.n_slices)} slices, darkness z {r.darkness_z}"
                          + (f"   [{', '.join(flags)}]" if flags else "") + feat
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
        save_calls(base, subj, df, calls)

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
not_rev = len(df) - len(calls)
if not_rev:
    print(f"WARNING: {not_rev} of {len(df)} candidates NOT reviewed - the score is a lower bound "
          f"(continue with: review_css.py {subj})")
cmd = [sys.executable, f"{base}/scripts/mark_css.py", subj, ",".join(map(str, accepted)) or "none",
       "--reviewed", ",".join(map(str, sorted(calls))) or "none"]
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
usage: mark_css.py SUBJECT 1,4,9 [--ich 3,7] [--reviewed 1-40 list]     (use 'none' when no cSS)
--reviewed: ids the reader actually looked at (review_css.py passes this). Default: all candidates.
--ich: candidates you judged to be siderosis/hemosiderin CONNECTED TO a lobar ICH.
       They are excluded from the multifocality score (standard convention) but reported."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, subprocess, pandas as pd
if len(sys.argv) < 3:
    sys.exit("usage: mark_css.py SUBJECT 1,4,9 [--ich 3,7]   (or: none)")
subj, arg = sys.argv[1], sys.argv[2]
ich_arg = sys.argv[sys.argv.index("--ich") + 1] if "--ich" in sys.argv else ""
rev_arg = sys.argv[sys.argv.index("--reviewed") + 1] if "--reviewed" in sys.argv else None
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
csv = f"{base}/review/{subj}_candidates.csv"
df = pd.read_csv(csv)
parse = lambda a: [] if a.lower() in ("", "none") else [int(x) for x in a.split(",") if x.strip()]
ids, ich = parse(arg), parse(ich_arg)
rev = parse(rev_arg) if rev_arg is not None else df.cand_id.astype(int).tolist()
valid = set(df.cand_id)
bad = [i for i in ids + ich + rev if i not in valid]
if bad: sys.exit(f"candidate ids not in list: {bad} (valid 1-{df.cand_id.max()})")
both = sorted(set(ids) & set(ich))
if both: sys.exit(f"ids marked both cSS and near-ICH: {both}")
df["accept"] = df.cand_id.isin(ids).astype(int)
df["near_ich"] = df.cand_id.isin(ich).astype(int)
df["reviewed"] = df.cand_id.isin(set(rev) | set(ids) | set(ich)).astype(int)
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
Candidates marked near_ich=1 (siderosis connected to a lobar ICH) are excluded and reported.
v4: SULCAL scoring when work/ID_a2009s_swispace.nii.gz exists (recon-all + prep_anat.sh):
  each accepted focus is assigned to the Destrieux sulci it lines (within 4 mm); two sulci are
  "immediately adjacent" if they touch or border the same gyrus. Per hemisphere (Charidimou):
  0 none; 1 = one sulcus or <=3 adjacent sulci; 2 = >=2 non-adjacent or >3 sulci.
  STRIVE-2 category: focal = 1-3 sulci, disseminated = >3 sulci.
  Without Destrieux labels the v3 Euclidean approximation is used (3 mm foci, 10 mm adjacency).
Infratentorial candidates (classical superficial siderosis pattern) are reported, not scored.
v4.3 ICH rule (Charidimou et al., Neurology 2017;89:2128): cSS "contiguous or potentially anatomically
connected with any lobar ICH" is not scored; cSS must be separated from any lobar ICH by >=3
unaffected sulci, or by >=2 (at multiple axial levels) if the haematoma has no superficial path along
the convexity. With Destrieux labels + an ICH mask, unaffected sulci between each accepted focus and
the ICH are counted on the sulcal adjacency graph: <2 -> excluded (reader-drawn ICH mask) or warned
(automatic mask); exactly 2 -> kept but flagged for the reader to check the 2-sulci conditions."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, json, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import load_lut

subj = sys.argv[1]
MERGE_MM, ADJ_MM, GROW_MM = 3.0, 10.0, 5.0
SULC_MM = 4.0                     # a focus lines a sulcus if within 4 mm of that sulcus' cortex
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

# ---- Destrieux sulcal labels (optional)
a2_p = f"{base}/work/{subj}_a2009s_swispace.nii.gz"
A2 = nib.load(a2_p).get_fdata().astype(np.int32) if os.path.exists(a2_p) else None
lut = load_lut()
is_sulc = lambda c: 11100 <= c < 12200 and ("_S_" in lut.get(c, "") or "Lat_Fis" in lut.get(c, ""))
is_gyr = lambda c: 11100 <= c < 12200 and not is_sulc(c)
if A2 is not None:
    SULC = np.array([c for c in np.unique(A2) if is_sulc(int(c))], np.int32)
    if not len(SULC):
        print("WARNING: Destrieux file has no sulcal labels known to the LUT - using Euclidean scoring")
        A2 = None

def crop(mask, mm):
    idx = np.argwhere(mask); pad = int(np.ceil(mm / min(vox))) + 1
    lo = np.maximum(idx.min(0) - pad, 0); hi = idx.max(0) + pad + 1
    return tuple(slice(l, h) for l, h in zip(lo, hi))

def sulci_of(mask, hemi):
    """Destrieux sulci lined by this focus: nearest sulcal-cortex label of each voxel within
    SULC_MM; labels holding >=15 % of those voxels (min 3) count."""
    sl = crop(mask, SULC_MM + 6); a = A2[sl]
    lo, hi = (11100, 11200) if hemi == "L" else (12100, 12200)
    sm = np.isin(a, SULC) & (a >= lo) & (a < hi)
    if not sm.any(): return []
    d, ind = ndi.distance_transform_edt(~sm, sampling=vox, return_indices=True)
    m = mask[sl]; near = m & (d <= SULC_MM)
    if not near.any(): near = m                      # gyral crown: take the nearest sulcus
    codes = a[tuple(i[near] for i in ind)]
    cnt = np.bincount(codes - lo, minlength=100)
    return [int(c + lo) for c in np.nonzero(cnt >= max(3, 0.15 * len(codes)))[0]]

_touch = {}
def touching(code):
    if code not in _touch:
        m = A2 == code; sl = crop(m, 3)
        dil = ndi.binary_dilation(m[sl], structure=np.ones((3, 3, 3), bool), iterations=2)
        _touch[code] = set(int(c) for c in np.unique(A2[sl][dil])) - {0, code}
    return _touch[code]

def adjacent(s, t):
    """immediately adjacent sulci: touch each other, or border the same gyrus"""
    ts, tt = touching(s), touching(t)
    return t in ts or s in tt or any(is_gyr(g) for g in ts & tt)

def components(nodes, edge):
    parent = {n: n for n in nodes}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i in nodes:
        for j in nodes:
            if i < j and edge(i, j): parent[find(i)] = find(j)
    return len({find(n) for n in nodes})

infra = df[df.infratentorial == 1].cand_id.astype(int).tolist() if "infratentorial" in df.columns else []
acc_all = df[df.accept == 1] if len(df) else df
acc = acc_all[~acc_all.cand_id.isin(infra)] if len(acc_all) else acc_all
acc_infra = acc_all[acc_all.cand_id.isin(infra)] if len(acc_all) else acc_all
# ---- ICH separation counted in sulci (Charidimou 2017)
ICH, ich_src = None, ""
for pth, src_ in ((f"{base}/work/{subj}_ich.nii.gz", "drawn"), (f"{base}/work/{subj}_ich_used.nii.gz", "auto")):
    if os.path.exists(pth):
        m_ = nib.load(pth).get_fdata() > 0
        ICH, ich_src = (m_, src_) if m_.any() else (None, src_)
        break
hemi_rng = lambda h: (11100, 11200) if h == "L" else (12100, 12200)

def ich_sulci():
    """Destrieux sulci next to the ICH (within 5 mm; else the nearest one)"""
    sl = crop(ICH, 12); a = A2[sl]; sm = np.isin(a, SULC)
    if not sm.any(): return set()
    d, ind = ndi.distance_transform_edt(~sm, sampling=vox, return_indices=True)
    m = ICH[sl]; near = m & (d <= 5.0)
    if not near.any(): near = m & (d <= d[m].min() + 1e-6)
    return set(int(c) for c in np.unique(a[tuple(i[near] for i in ind)]))

def sulcal_hops(start, h):
    """breadth-first distance (in sulcal steps) on the adjacency graph of one hemisphere"""
    lo, hi = hemi_rng(h); nodes = [int(c) for c in SULC if lo <= c < hi]
    dist = {c: 0 for c in start if lo <= c < hi}; front = list(dist)
    while front:
        nxt = []
        for s_ in front:
            for t in nodes:
                if t not in dist and adjacent(s_, t): dist[t] = dist[s_] + 1; nxt.append(t)
        front = nxt
    return dist

cand_sulci, ich_excl, ich_check, ich_between = {}, [], [], {}
if A2 is not None and len(acc):
    for _, r in acc.iterrows():
        cand_sulci[int(r.cand_id)] = sulci_of(lab == int(r.cand_id), r.hemi)
    if ICH is not None:
        IS = ich_sulci(); hops = {h: sulcal_hops(IS, h) for h in ("L", "R")}
        for _, r in acc.iterrows():
            c = int(r.cand_id); dd = [hops[r.hemi][x] for x in cand_sulci[c] if x in hops[r.hemi]]
            if not dd: continue                       # not connected to the ICH's sulci (e.g. other hemisphere)
            between = min(dd) - 1                     # unaffected sulci in between (-1 = same sulcus)
            ich_between[c] = between
            if between < 2: ich_excl.append(c)
            elif between == 2: ich_check.append(c)
        if ich_src == "drawn" and ich_excl:
            acc = acc[~acc.cand_id.isin(ich_excl)]

method = "sulcal (Destrieux)" if A2 is not None else "euclidean (approximation)"
result = {"subject": subj, "scoring_method": method}; total = total_foci = total_sulci = 0
sulci_names = []
for h in ["L", "R"]:
    ids = acc[acc.hemi == h].cand_id.astype(int).tolist() if len(acc) else []
    masks = {i: lab == i for i in ids}
    D = {(i, j): min_dist(masks[i], masks[j]) for i in ids for j in ids if i < j}
    foci = groups(ids, D, MERGE_MM); clusters = groups(ids, D, ADJ_MM)
    fpc = [sum(1 for f in foci if f[0] in c) for c in clusters]
    if A2 is not None:
        S = sorted(set(c for i in ids for c in cand_sulci.get(i) or sulci_of(masks[i], h)))
        ncomp = components(S, adjacent) if S else 0
        score = 0 if not ids else (1 if len(S) <= 3 and ncomp <= 1 else 2)
        result[f"{h}_sulci"] = len(S); total_sulci += len(S)
        # internally calibrated extent (van Harten et al. 2023 suggest % of sulci affected as less
        # sequence-dependent than volume): affected / all Destrieux sulci of the hemisphere
        n_all = int(sum(1 for c in SULC if hemi_rng(h)[0] <= c < hemi_rng(h)[1]))
        result[f"{h}_sulci_pct"] = round(100.0 * len(S) / max(n_all, 1), 1)
        sulci_names += [lut.get(c, str(c)) for c in S]
    else:
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
n_units = total_sulci if A2 is not None else total_foci     # STRIVE-2 counts sulci
result["category"] = "absent" if n_units == 0 else ("focal" if n_units <= 3 else "disseminated")
if A2 is not None:
    result["n_sulci_total"] = total_sulci
    result["sulci"] = sulci_names
result["infratentorial_accepted"] = int(len(acc_infra))
if A2 is not None and ICH is not None:
    result["ich_mask"] = ich_src
    result["ich_sulcal_too_close"] = ",".join(map(str, sorted(ich_excl)))
    result["ich_sulcal_check"] = ",".join(map(str, sorted(ich_check)))
# sequence (SWI vs T2*-GRE scores are not interchangeable: SWI rates higher)
seq = "unknown"
import glob as _glob
for jp in sorted(_glob.glob(f"{base}/raw/{subj}/SWI/*.json")):
    try:
        it = " ".join(json.load(open(jp)).get("ImageType", [])).upper()
        seq = "SWI" if "SWI" in it else ("T2*-GRE" if it else seq)
        break
    except Exception:
        pass
result["sequence"] = seq
result["candidate_volume_mm3"] = round(float(acc.volume_mm3.sum()) if len(acc) else 0.0, 1)
result["grown_volume_mm3"] = round(grown_vol, 1)
result["regions"] = sorted(acc.region.unique().tolist()) if len(acc) else []
ich = df[df.near_ich == 1] if "near_ich" in df.columns and len(df) else df.iloc[0:0]
result["n_candidates"] = int(len(df))
result["n_unreviewed"] = int((df.reviewed == 0).sum()) if "reviewed" in df.columns and len(df) else 0
result["near_ich_candidates"] = int(len(ich))
result["near_ich_volume_mm3"] = round(float(ich.volume_mm3.sum()) if len(ich) else 0.0, 1)

json.dump(result, open(f"{base}/review/{subj}_score.json", "w"), indent=2)
summ = f"{base}/review/css_scores.csv"
row = pd.DataFrame([{k: v for k, v in result.items() if k not in ("regions", "sulci")}])
if os.path.exists(summ):
    old = pd.read_csv(summ); old = old[old.subject != subj]
    row = pd.concat([old, row], ignore_index=True)
row.to_csv(summ, index=False)

print(f"\n{subj}  cSS multifocality score: {total}/4   ({result['category']})   [{method}; {seq}]")
for h, nm in (("L", "left "), ("R", "right")):
    print(f"  {nm}: score {result[f'{h}_score']}  foci {result[f'{h}_foci']}  clusters {result[f'{h}_clusters']}"
          + (f"  sulci {result[f'{h}_sulci']} ({result[f'{h}_sulci_pct']}% of sulci)" if A2 is not None else ""))
if A2 is not None and sulci_names:
    print(f"  sulci: {', '.join(sulci_names)}")
print(f"  volume: candidates {result['candidate_volume_mm3']} mm3, grown (full extent) {result['grown_volume_mm3']} mm3")
print(f"  regions: {', '.join(result['regions'])}")
if result["n_unreviewed"]:
    print(f"  WARNING: {result['n_unreviewed']} of {result['n_candidates']} candidates were not reviewed "
          f"- score is a LOWER BOUND")
if result["infratentorial_accepted"]:
    print(f"  infratentorial siderosis (NOT in cSS score - consider classical superficial siderosis): "
          f"{result['infratentorial_accepted']} candidates")
if ich_excl:
    print(f"  ICH rule (<3 unaffected sulci to the lobar ICH, Charidimou 2017): candidates {sorted(ich_excl)} "
          + ("EXCLUDED from the score" if ich_src == "drawn" else
             "would be excluded - automatic ICH mask: confirm the ICH (draw work/" + subj + "_ich.nii.gz) "
             "or call them 'i' in the review"))
if ich_check:
    print(f"  ICH rule: candidates {sorted(ich_check)} are exactly 2 sulci from the ICH - kept; exclude ('i') if "
          f"they are not 2 sulci away at multiple axial levels or the haematoma reaches the surface")
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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, load_calls
base = _base()
rows = []
for f in glob.glob(f"{base}/review/*_calls.csv"):
    s = os.path.basename(f).replace("_calls.csv", "")
    d = pd.read_csv(f"{base}/review/{s}_candidates.csv")
    calls, warn = load_calls(base, s, d)          # v4: matched by fingerprint, not by rank
    if warn: print("WARNING:", warn)
    if not calls: continue
    m = d.merge(pd.DataFrame(list(calls.items()), columns=["cand_id", "call"]), on="cand_id")
    m["subject"] = s; rows.append(m)
if not rows: sys.exit("no reviewed cases yet (run review_css.py first)")
D = pd.concat(rows, ignore_index=True)
D = D[D.call.isin(["cSS", "Vein", "Normal", "Artifact"])]
print(f"{D.subject.nunique()} subject(s), {len(D)} decided candidates: {D.call.value_counts().to_dict()}\n")
def auc(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float); x, y = x[~np.isnan(x)], y[~np.isnan(y)]
    if not len(x) or not len(y): return float("nan")
    return float((x[:, None] > y[None, :]).mean() + 0.5 * (x[:, None] == y[None, :]).mean())
feats = ["darkness_z", "volume_mm3", "n_slices", "elongation", "surface_contact", "cortex_frac",
         "cortex_dist_mm", "branch_per10mm", "surface_gradient",
         "pial_dist_mm", "bank_frac", "tube_ratio", "surface_alignment", "tram_frac", "vein_tree_mm",
         "mirror_dark_frac", "flair_csf_z", "score_v3", "score_v4", "score"]
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
MIN_COVER = 0.30   # v4: a lesion counts as found only if candidates cover >=30 % of it (1-voxel
                   # tolerance for blooming); v3 accepted any touching candidate, e.g. a vein crossing it
n = int(truth.max()); best, lesion_c = [], set()
for i in range(1, n + 1):
    t = truth == i; nt = int(t.sum())
    touch = sorted(set(np.unique(cand[ndi.binary_dilation(t, iterations=1)])) - {0})
    cov = {c: (ndi.binary_dilation(cand == c, iterations=1) & t).sum() / nt for c in touch}
    ids = [c for c in touch if cov[c] >= 0.05]
    total = float((ndi.binary_dilation(np.isin(cand, ids), iterations=1) & t).sum() / nt) if ids else 0.0
    found = total >= MIN_COVER
    if found: best.append(ids[0]); lesion_c |= set(ids)
    print(f"  lesion {i}: " + (f"FOUND  (candidate #{ids[0]}, coverage {total:.0%})" if found
                               else f"MISSED (coverage {total:.0%})"))
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
              "surface_gradient", "surface_contact", "elongation", "pial_dist_mm", "bank_frac",
              "tube_ratio", "surface_alignment", "tram_frac", "vein_tree_mm", "mirror_dark_frac",
              "score_v3", "score_v4", "score"]:
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

## tests/phantom_v4.py

```python
"""Phantom brain with the structures the cSS definition separates (v4 tests).
usage: phantom_v4.py NAME SEED        (writes into $CSS_BASE; appends test labels to the LUT)

Geometry (mm, 0.8 x 0.8 x 1.2 voxels like a clinical SWI):
  ellipsoid hemispheres, 2.5 mm cortex, 4.5 mm extra-axial CSF (SynthStrip keeps CSF), midline
  fissure, 16 radial sulci (2 mm CSF, 14 mm deep).
Truth  work/NAME_truth.nii.gz  (cSS, ids 1-6):
  1-3 tram-track: subpial band on BOTH banks of a sulcus        (the classic sign)
  4-5 convexity: subpial band following the gyral crown          (curvilinear sheet)
  6   single bank of a sulcus
Mimics work/NAME_veins.nii.gz (ids 1-6):
  1-4 tubular veins in the middle of the sulcal CSF (1,2 and 3,4 are left/right mirror pairs)
  5-6 cortical surface veins running tangentially in the CSF just outside the pial surface
Also: a lobar ICH (dark sphere, 8 mm) in left white matter -> work/NAME_ich_truth.nii.gz
Labels: synthseg/NAME_seg.nii.gz (aseg/DK-like) and work/NAME_a2009s_swispace.nii.gz
(Destrieux-like: S_k = banks of sulcus k, G_k = gyrus between sulci k and k+1)."""
import os, sys, numpy as np, nibabel as nib
from scipy import ndimage as ndi

B = os.environ["CSS_BASE"]
name, sd = sys.argv[1], int(sys.argv[2]); rng = np.random.default_rng(sd)
vox = (0.8, 0.8, 1.2); sh = (150, 165, 90)
aff = np.diag(list(vox) + [1.0])
x, y, z = np.meshgrid(*[(np.arange(n) - n / 2) * v for n, v in zip(sh, vox)], indexing="ij")
R0 = 52.0
r = np.sqrt((x / 52) ** 2 + (y / 58) ** 2 + (z / 44) ** 2)
depth = (1 - r) * R0                        # ~mm below the pial surface (negative = outside)
inside = r < 1                              # pial surface
mask = depth > -4.5                         # brain mask incl. 4.5 mm extra-axial CSF (atrophic CAA)
left = x < 0
rho = np.sqrt(x ** 2 + y ** 2); ang = np.arctan2(y, x)

seg = np.zeros(sh, np.int32)
seg[mask] = 24
seg[inside] = np.where(left, 2, 41)[inside]
ctx = np.where(left, 1000, 2000) + np.where(y > 0, 28, 22)
cort = inside & (depth < 2.5)
A0 = np.arange(-np.pi, np.pi, np.pi / 8) + np.pi / 16          # 16 sulci, none on the midline
sul_any = np.zeros(sh, bool); bank_k = np.full(sh, -1, np.int16)
dl = np.zeros((len(A0),) + sh, np.float32)
for k, a0 in enumerate(A0):
    d_line = np.abs(np.sin(ang - a0)) * rho
    front = np.cos(ang - a0) > 0
    dl[k] = np.where(front, d_line, 99)
    sul = inside & (dl[k] < 1.0) & (depth < 14)
    bank = inside & (dl[k] < 3.5) & (depth < 16.5) & ~sul
    cort |= bank; sul_any |= sul
    bank_k[bank & (bank_k < 0)] = k
fiss = mask & (np.abs(x) < 1.0)
cort |= inside & (np.abs(x) < 3.5) & ~fiss
cort &= ~sul_any & ~fiss
seg[cort] = ctx[cort]; seg[sul_any | fiss] = 24
seg[~mask] = 0

# Destrieux-like labels on the cortex
sector = (np.floor((ang - A0[0]) / (np.pi / 8)).astype(int)) % len(A0)
a2 = np.zeros(sh, np.int32)
hemi_off = np.where(left, 11100, 12100)
a2[cort] = (hemi_off + 1 + sector)[cort]                         # G_k
bk = cort & (bank_k >= 0)
a2[bk] = (hemi_off + 50 + bank_k)[bk]                            # S_k

I = np.zeros(sh, np.float32)
I[np.isin(seg, [2, 41])] = 180; I[seg >= 1000] = 210; I[seg == 24] = 250
I[mask] += rng.normal(0, 10, int(mask.sum()))

def darken(m, f=0.6, blur=0.6):
    """hemosiderin / deoxy-Hb: multiplicative darkening with blooming (gaussian edge)"""
    w = np.clip(ndi.gaussian_filter(m.astype(np.float32), [blur / v for v in vox]) * 1.6, 0, 1)
    I[:] = I * (1 - f * w)

truth = np.zeros(sh, np.int16); veins = np.zeros(sh, np.int16)
zc = lambda lo, hi: (z > lo) & (z < hi)
# 1-3 tram-track: subpial 1.2 mm band on both banks (sulci 2, 6 left side; 11 right side)
for i, (k, z0) in enumerate([(5, -6), (7, 4), (13, -10)], 1):
    m = (dl[k] >= 1.0) & (dl[k] < 2.2) & (depth > 3) & (depth < 12) & zc(z0, z0 + 16) & cort
    truth[m] = i
# 4-5 convexity crown: outer 1 mm of cortex over a 14 deg arc, 18 mm high
for i, (a, z0) in enumerate([(A0[4] + np.pi / 16, -2), (A0[12] + np.pi / 16, 6)], 4):
    m = cort & (depth >= 0) & (depth < 1.2) & (np.abs(np.angle(np.exp(1j * (ang - a)))) < np.radians(7)) & zc(z0, z0 + 18)
    truth[m] = i
# 6 single bank
k = 9; s = np.sign(np.sin(ang - A0[k]))
m = (dl[k] >= 1.0) & (dl[k] < 2.2) & (s > 0) & (depth > 3) & (depth < 12) & zc(-4, 12) & cort
truth[m] = 6
# veins 1-4: tubes (r 0.7 mm) in the sulcal CSF at 5 mm depth, running 30 mm along z;
# 1/2 and 3/4 are left-right mirror images (normal anatomy is roughly symmetric)
# (sulcus k mirrors sulcus 7-k mod 16; sulci 5, 7, 9, 13 carry cSS)
for i, k in enumerate([1, 6, 11, 12], 1):
    c_xy = np.array([np.cos(A0[k]), np.sin(A0[k])])
    # point on the sulcal plane where depth ~5 mm (search along the ray)
    t = np.linspace(20, 60, 400); px, py = c_xy[0] * t, c_xy[1] * t
    rr = np.sqrt((px / 52) ** 2 + (py / 58) ** 2); p = np.argmin(np.abs((1 - rr) * R0 - 5))
    m = (np.sqrt((x - px[p]) ** 2 + (y - py[p]) ** 2) < 0.7) & zc(-14, 16) & mask
    veins[m] = i
# veins 5-6: cortical surface veins 1.2 mm outside the pial surface, tangential, 40 deg arc
for i, (a, z0) in enumerate([(np.radians(60), 8), (np.radians(-120), -6)], 5):
    m = (np.abs(depth + 1.2) < 0.7) & (np.abs(z - z0) < 0.7) & \
        (np.abs(np.angle(np.exp(1j * (ang - a)))) < np.radians(20)) & mask
    veins[m] = i
ich = (np.sqrt((x + 22) ** 2 + (y - 10) ** 2 + (z - 4) ** 2) < 8) & inside
darken(truth > 0, 0.6); darken(veins > 0, 0.6, blur=0.4); darken(ich, 0.7, blur=0.8)
I[~mask] = 0; I[mask] = np.maximum(I[mask], 1)

os.makedirs(f"{B}/data", exist_ok=True); os.makedirs(f"{B}/work", exist_ok=True)
os.makedirs(f"{B}/synthseg", exist_ok=True)
nib.save(nib.Nifti1Image(I, aff), f"{B}/data/{name}_swi.nii")
nib.save(nib.Nifti1Image(seg, aff), f"{B}/synthseg/{name}_seg.nii.gz")
nib.save(nib.Nifti1Image(a2, aff), f"{B}/work/{name}_a2009s_swispace.nii.gz")
nib.save(nib.Nifti1Image(truth, aff), f"{B}/work/{name}_truth.nii.gz")
nib.save(nib.Nifti1Image(veins, aff), f"{B}/work/{name}_veins.nii.gz")
nib.save(nib.Nifti1Image(ich.astype(np.uint8), aff), f"{B}/work/{name}_ich_truth.nii.gz")
lutp = os.path.join(os.environ.get("FREESURFER_HOME", ""), "FreeSurferColorLUT.txt")
if os.path.exists(lutp) and "ctx_lh_S_k00" not in open(lutp).read():
    with open(lutp, "a") as f:
        for h, off in (("lh", 11100), ("rh", 12100)):
            for k in range(len(A0)):
                f.write(f"{off + 1 + k} ctx_{h}_G_k{k:02d} 100 100 100 0\n")
                f.write(f"{off + 50 + k} ctx_{h}_S_k{k:02d} 200 200 200 0\n")
print(f"{name}: phantom v4  cSS lesions {int(truth.max())}  veins {int(veins.max())}  "
      f"ICH {int(ich.sum() * np.prod(vox))} mm3")
```

## tests/check_v4.py

```python
"""Check v4 detector output on the phantom_v4 truth.   usage: check_v4.py NAME
Prints each candidate's truth class and features, sensitivity, and how well score_v3 / score_v4
and each definition feature separate true cSS from mimics. Exits 1 if a sanity check fails."""
import os, sys, numpy as np, nibabel as nib, pandas as pd
B = os.environ["CSS_BASE"]; name = sys.argv[1]
ld = lambda p: nib.load(p).get_fdata().astype(int)
truth, veins = ld(f"{B}/work/{name}_truth.nii.gz"), ld(f"{B}/work/{name}_veins.nii.gz")
cand = ld(f"{B}/review/{name}_candidates.nii.gz")
df = pd.read_csv(f"{B}/review/{name}_candidates.csv")
cls = []
for c in df.cand_id:
    m = cand == c; t = np.bincount(truth[m], minlength=7)[1:]; v = np.bincount(veins[m], minlength=7)[1:]
    if t.sum() >= v.sum() and t.sum() > 0.2 * m.sum(): cls.append(f"cSS{t.argmax() + 1}")
    elif v.sum() > 0: cls.append(f"vein{v.argmax() + 1}")
    else: cls.append("other")
df["truth"] = cls
F = ["darkness_z", "pial_dist_mm", "bank_frac", "tube_ratio", "surface_alignment", "tram_frac",
     "vein_tree_mm", "mirror_dark_frac", "score_v3", "score_v4"]
print(df[["cand_id", "truth"] + F].to_string(index=False))
found = sorted({int(c[3:]) for c in cls if c.startswith("cSS")})
n = int(truth.max()); print(f"\nsensitivity: {len(found)}/{n}  found lesions {found}")
pos = df[df.truth.str.startswith("cSS")]; neg = df[~df.truth.str.startswith("cSS")]
def auc(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float); x, y = x[~np.isnan(x)], y[~np.isnan(y)]
    if not len(x) or not len(y): return float("nan")
    return float((x[:, None] > y[None, :]).mean() + 0.5 * (x[:, None] == y[None, :]).mean())
print("AUC cSS vs mimics (>0.5 = higher in cSS):  " +
      "  ".join(f"{f} {auc(pos[f], neg[f]):.2f}" for f in F))
ich = nib.load(f"{B}/work/{name}_ich_used.nii.gz").get_fdata() > 0
icht = nib.load(f"{B}/work/{name}_ich_truth.nii.gz").get_fdata() > 0
ich_dice = 2 * (ich & icht).sum() / max(ich.sum() + icht.sum(), 1)
print(f"automatic ICH mask Dice vs truth: {ich_dice:.2f}")
ok = len(found) >= n - 1 and auc(pos.score_v4, neg.score_v4) >= 0.8 and ich_dice > 0.5
print("RESULT v4 " + ("OK" if ok else "FAIL") +
      f" found={len(found)} total={n} auc_v3={auc(pos.score_v3, neg.score_v3):.2f} "
      f"auc_v4={auc(pos.score_v4, neg.score_v4):.2f} ich_dice={ich_dice:.2f}")
sys.exit(0 if ok else 1)
```

## tests/run_tests.sh

```bash
#!/bin/bash
# Self-test with a synthetic phantom brain - no patient data, no FreeSurfer needed.
# usage: bash tests/run_tests.sh      (from the package root, inside the 'css' conda env)
# Shell flow of run_css.sh / prep_anat.sh with stubbed FreeSurfer: bash tests/test_shell.sh
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
echo "== calls survive a detector re-run (fingerprint matching)"
cp $T/review/PH1S_candidates.csv $T/old_cand.csv
python $S/detect_css.py PH1S -2.3 80 3 > /dev/null
python - "$T" <<'PYEOF'
import sys, pandas as pd; sys.path.insert(0, sys.argv[1] + "/scripts")
from css_common import load_calls, FP
b = sys.argv[1]; old = pd.read_csv(b + "/old_cand.csv"); new = pd.read_csv(b + "/review/PH1S_candidates.csv")
raw = pd.read_csv(b + "/review/PH1S_calls.csv"); calls, warn = load_calls(b, "PH1S", new)
for cid, call in calls.items():          # every re-attached call must sit on the same lesion as before
    r = raw[raw.call == call]; n = new.set_index("cand_id").loc[cid]
    assert (((r[FP[:3]] - n[FP[:3]].astype(float)) ** 2).sum(axis=1) ** .5 <= 2).any(), (cid, call)
print(f"   {len(calls)}/{len(raw)} calls re-attached by fingerprint" + (f"; {warn}" if warn else ""))
PYEOF
python $S/detect_css.py PH1S > /dev/null
echo "== feature report";      python $S/feature_report.py | head -3
echo "== stress test (small)"; python $S/stress_test.py --hosts PH1,PH2 --depths 0.6 --seeds 1 --tag test | grep -A2 SUMMARY
echo "== v4 phantom: tram-track / convexity cSS vs tubular & surface veins, ICH, sulcal scoring"
python $ROOT/tests/phantom_v4.py PH3 1; python $S/align_seg.py PH3 > /dev/null
python $S/detect_css.py PH3 | head -1
python $ROOT/tests/check_v4.py PH3 | tail -3
python $S/score_css.py PH3 --truth > $T/score_v4.log
grep -q "3/4" $T/score_v4.log && grep -q "sulcal (Destrieux)" $T/score_v4.log || { cat $T/score_v4.log; echo "FAIL: sulcal score on PH3 should be 3/4 (L1 + R2)"; exit 1; }
grep -A3 "cSS multifocality" $T/score_v4.log
echo "== ICH rule in sulci (Charidimou 2017): reader-drawn ICH -> left foci next to it excluded"
cp $T/work/PH3_ich_truth.nii.gz $T/work/PH3_ich.nii.gz
python $S/score_css.py PH3 --truth > $T/score_ich.log; rm $T/work/PH3_ich.nii.gz
grep -q "2/4" $T/score_ich.log && grep -q "EXCLUDED" $T/score_ich.log || { cat $T/score_ich.log; echo "FAIL: ICH sulcal rule"; exit 1; }
grep "ICH rule" $T/score_ich.log
echo; echo "ALL TESTS PASSED  (temp dir $T)"
```

## tests/test_shell.sh

```bash
#!/bin/bash
# Shell-flow test for run_css.sh + prep_anat.sh with STUBBED FreeSurfer tools and python
# (checks argument parsing, which tools run in which order, and the re-use / staleness logic).
# usage: bash tests/test_shell.sh        Uses real zsh if installed, else bash as a stand-in.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
T=$(mktemp -d); mkdir -p $T/bin $T/base/scripts; cp $ROOT/scripts/*.py $ROOT/scripts/*.sh $T/base/scripts/
cat > $T/bin/fsstub <<'EOS'
#!/bin/bash
n=$(basename $0); echo "$n $*" >> $STUBLOG
out=""; in=""; targ=""
while [ $# -gt 0 ]; do case "$1" in
  -o|--o) out=$2; shift 2;; -i|--i|--mov) [ -z "$in" ] && in=$2; shift 2;; --targ) targ=$2; shift 2;;
  --reg) touch $2; shift 2;;
  -s) [ "$n" = recon-all ] && mkdir -p $SUBJECTS_DIR/$2/mri && touch $SUBJECTS_DIR/$2/mri/aparc+aseg.mgz $SUBJECTS_DIR/$2/mri/aparc.a2009s+aseg.mgz; shift 2;;
  *) [ "$n" = mri_convert ] && { [ -z "$in" ] && in=$1 || out=$1; }; shift;; esac; done
[ -n "$out" ] && [ -n "$in" ] && cp "$in" "$out" 2>/dev/null || true
EOS
chmod +x $T/bin/fsstub
for t in mri_convert mri_synthstrip mri_synthseg mri_coreg mri_vol2vol bbregister recon-all freeview; do ln -s fsstub $T/bin/$t; done
printf '#!/bin/bash\necho "python $(basename $1) ${*:2}" >> $STUBLOG\n' > $T/bin/python; chmod +x $T/bin/python
command -v zsh >/dev/null || { printf '#!/bin/bash\nexec bash "$@"\n' > $T/bin/zsh; chmod +x $T/bin/zsh; echo "(zsh not installed: using bash as a stand-in)"; }
echo swi > $T/swi.nii; echo t1 > $T/t1.nii; echo fl > $T/fl.nii
export PATH=$T/bin:$PATH STUBLOG=$T/log CSS_BASE=$T/base NOVIEW=1 SUBJECTS_DIR=/somewhere/else  # must be ignored
SH="zsh"
fail() { echo "FAIL: $1"; cat $T/log; exit 1; }
seen() { grep -q "$1" $T/log || fail "expected call: $1"; }
notseen() { ! grep -q "$1" $T/log || fail "unexpected call: $1"; }

echo "== 1. SWI only"
: > $T/log; $SH $T/base/scripts/run_css.sh P1 $T/swi.nii > /dev/null
seen "mri_synthstrip -i $T/base/work/P1_raw"; seen "mri_synthseg --i $T/base/data/P1_swi.nii"
seen "python align_seg.py P1"; seen "python detect_css.py P1 -2.5 85 3"; notseen mri_coreg
echo "== 2. re-run: SynthSeg re-used"
: > $T/log; $SH $T/base/scripts/run_css.sh P1 > /dev/null; notseen mri_synthseg; seen "detect_css.py P1"
echo "== 3. SWI re-imported -> stale SynthSeg is redone"
sleep 1; : > $T/log; $SH $T/base/scripts/run_css.sh P1 $T/swi.nii > /dev/null; seen "mri_synthseg --i $T/base/data/P1_swi.nii"
echo "== 4. SWI + T1 + FLAIR (SynthSeg on T1, mri_coreg)"
: > $T/log; $SH $T/base/scripts/run_css.sh P2 $T/swi.nii --t1 $T/t1.nii --flair $T/fl.nii > /dev/null
seen "mri_synthseg --i $T/base/data/P2_t1.nii.gz"; seen "mri_coreg --mov $T/base/data/P2_swi.nii --ref $T/base/work/P2_t1_brain"
seen "mri_vol2vol --mov $T/base/data/P2_swi.nii --targ $T/base/work/P2_t1synthseg.nii.gz --lta $T/base/work/P2_swi2t1.lta --inv --interp nearest"
seen "mri_coreg --mov $T/base/work/P2_flair_brain.nii.gz --ref $T/base/data/P2_swi.nii"
seen "mri_vol2vol --mov $T/base/work/P2_flair_brain.nii.gz --targ $T/base/data/P2_swi.nii"
notseen "mri_synthseg --i $T/base/data/P2_swi.nii"; notseen recon-all; notseen bbregister
echo "== 5. --recon (recon-all, bbregister, Destrieux)"
: > $T/log; $SH $T/base/scripts/run_css.sh P3 $T/swi.nii --t1 $T/t1.nii --recon > /dev/null
seen "recon-all -s P3"; seen "bbregister --s P3 --mov $T/base/data/P3_swi.nii --reg $T/base/work/P3_swi2t1.lta --t2 --init-coreg"
seen "aparc.a2009s+aseg.mgz --lta $T/base/work/P3_swi2t1.lta --inv --interp nearest --o $T/base/work/P3_a2009s_swispace"
echo "== 6. recon-all already done -> re-used without --recon"
: > $T/log; $SH $T/base/scripts/prep_anat.sh P3 --t1 $T/t1.nii > /dev/null; notseen "recon-all"; seen bbregister
echo "== 7. batch run_all.sh reaches detection for every subject"
: > $T/log; $SH $T/base/scripts/run_all.sh P1 P2 > /dev/null || true   # stub python prints no summary
[ $(grep -c "detect_css.py" $T/log) -eq 2 ] || fail "run_all.sh did not run detection for both subjects"
echo "SHELL TESTS PASSED"
```
