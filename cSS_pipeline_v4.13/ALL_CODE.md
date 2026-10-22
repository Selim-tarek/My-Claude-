# All code — cSS pipeline v4.19

Generated from the files below; the files are authoritative.

## scripts/run_css.sh

```bash
#!/bin/zsh
# One patient, start to finish (up to review).
# usage: run_css.sh SUBJECT_ID [path/to/swi.nii(.gz)] [--nostrip] [--t1 T1.nii] [--flair FLAIR.nii] [--recon]
#                   [--phase PHASE.nii]
#   - with a SWI path: imports the scan; skull-strips it with SynthStrip unless --nostrip
#     (use --nostrip only for scans that are already skull-stripped)
#   - SWI = the processed SWI (or magnitude) series, NOT the minIP and NOT the phase map
#   - --t1 / --flair / --recon: v4 anatomy from T1 (and FLAIR), see prep_anat.sh
#     (T1 labels are re-used on later runs; SynthSeg on the SWI is then skipped)
#   - --phase: the SWI PHASE series of the same acquisition (for reading only: the reviewer shows the
#     filtered phase, where calcium has the opposite sign to veins/blood products)
#   - set NOVIEW=1 to skip opening freeview (batch use)
set -e
S=$1; shift 2>/dev/null || true
SRC=""; STRIP=1; ANAT=(); PHASE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --nostrip) STRIP=0; shift;;
    --t1|--flair) ANAT+=("$1" "$2"); shift 2;;
    --recon) ANAT+=("$1"); shift;;
    --phase) PHASE=$2; shift 2;;
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
    mri_convert "$SRC" $B/work/${S}_raw.nii.gz > /dev/null || { echo "ERROR: mri_convert could not read $SRC"; exit 1; }
    # v4.18.3: errors were silent under set -e (LEMON P203 stopped here with no message); keep SynthStrip's output
    mri_synthstrip -i $B/work/${S}_raw.nii.gz -o $B/data/${S}_swi.nii \
      -m $B/work/${S}_brainmask.nii.gz > $B/work/${S}_synthstrip.log 2>&1 \
      || { echo "ERROR: mri_synthstrip failed on the SWI (exit $?) - last lines of $B/work/${S}_synthstrip.log:"; tail -5 $B/work/${S}_synthstrip.log; exit 1; }
  else
    echo "[0/4] importing (no skull-strip) $SRC"
    mri_convert "$SRC" $B/data/${S}_swi.nii > /dev/null
  fi
fi
if [ ! -f $B/data/${S}_swi.nii ]; then
  echo "ERROR: $B/data/${S}_swi.nii not found - give the scan path as 2nd argument"; exit 1
fi
if [ -n "$PHASE" ]; then
  echo "[0a] importing SWI phase (for the reviewer)"
  mri_convert "$PHASE" $B/data/${S}_phase.nii.gz > /dev/null
fi
if [ ${#ANAT[@]} -gt 0 ]; then
  echo "[0b] T1/FLAIR anatomy -> SWI grid"
  zsh $SC/prep_anat.sh $S "${ANAT[@]}" || { echo "ERROR: prep_anat.sh failed - stopping (old labels would be used)"; exit 1; }
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
    --parc --robust --threads ${CSS_THREADS:-8} > /dev/null
else
  echo "      SynthSeg on the SWI (already done, skipping)"
fi
echo "[2/4] aligning segmentation to SWI"
python $SC/align_seg.py $S
echo "[3/4] detecting cSS candidates"
python $SC/detect_css.py $S -2.5 85 3
# v4.19: pial-surface features (incl. curve-vs-patch width) whenever recon-all surfaces exist; recorded only,
# and a failure here (e.g. the surface/label registration check) never stops the run
SD=${CSS_SUBJECTS_DIR:-$B/subjects}
if [ -f $SD/$S/surf/lh.pial ] && [ -f $SD/$S/mri/orig.mgz ] && [ -f $B/work/${S}_swi2t1.lta ]; then
  echo "[3b] pial-surface features (recon-all surfaces; recorded only)"
  python $SC/surface_features.py $S || echo "      WARNING: surface_features.py failed - detection results are unaffected"
fi
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
# v4.17.2: stop at the first failed step (a typo in --t1 used to fall through and leave the old labels in place)
die() { echo "ERROR: $1 - nothing after this step was run; the previous labels are unchanged"; exit 1; }
[ -f "$T1" ] || die "T1 file not found: $T1"
B=${CSS_BASE:-$HOME/css_project}; W=$B/work
# always the project folder (FreeSurfer's setup script points SUBJECTS_DIR at its own install dir);
# set CSS_SUBJECTS_DIR to use a different folder
export SUBJECTS_DIR=${CSS_SUBJECTS_DIR:-$B/subjects}
SWI=$B/data/${S}_swi.nii
[ -f $SWI ] || { echo "ERROR: $SWI not found - run run_css.sh $S <swi> first"; exit 1; }
mkdir -p $W $SUBJECTS_DIR

echo "[a] T1: import + skull-strip (SynthStrip)"
# via a temporary file, so --t1 data/ID_t1.nii.gz (the copy imported earlier) also works
mri_convert "$T1" $B/data/${S}_t1_tmp.nii.gz > /dev/null || die "mri_convert could not read $T1"
mv $B/data/${S}_t1_tmp.nii.gz $B/data/${S}_t1.nii.gz
mri_synthstrip -i $B/data/${S}_t1.nii.gz -o $W/${S}_t1_brain.nii.gz > /dev/null || die "mri_synthstrip failed"

A2=""
if [ $RECON -eq 1 ] || [ -f $SUBJECTS_DIR/$S/mri/aparc.a2009s+aseg.mgz ]; then
  if [ ! -f $SUBJECTS_DIR/$S/mri/aparc.a2009s+aseg.mgz ]; then
    # v4.18: an unfinished subject folder (earlier failed run) blocks recon-all -i; keep it, renamed
    if [ -d $SUBJECTS_DIR/$S ]; then
      mv $SUBJECTS_DIR/$S $SUBJECTS_DIR/${S}_failed_$(date +%Y%m%d_%H%M%S)
      echo "    moved the unfinished $SUBJECTS_DIR/$S aside (${S}_failed_<date>)"
    fi
    echo "[b] recon-all -all (1-3 h; log: $W/${S}_recon.log)"
    recon-all -s $S -i $B/data/${S}_t1.nii.gz -all -threads ${CSS_THREADS:-8} > $W/${S}_recon.log 2>&1 \
      || die "recon-all failed (see $W/${S}_recon.log)"
  else
    echo "[b] recon-all already finished for $S - re-using it"
  fi
  LAB=$SUBJECTS_DIR/$S/mri/aparc+aseg.mgz
  A2=$SUBJECTS_DIR/$S/mri/aparc.a2009s+aseg.mgz
  echo "[c] register SWI -> T1 (bbregister, boundary-based, T2*-like contrast)"
  bbregister --s $S --mov $SWI --reg $W/${S}_swi2t1.lta --t2 --init-coreg > $W/${S}_bbreg.log 2>&1 \
    || die "bbregister failed (see $W/${S}_bbreg.log)"
  # FreeSurfer 8 bbregister names its side files after an internal .dat (ID_swi2t1.dat.mincost), older ones .lta.mincost
  MC=$(ls $W/${S}_swi2t1.lta.mincost $W/${S}_swi2t1.dat.mincost 2>/dev/null | head -1)
  REGM=bbregister; COST=$( [ -n "$MC" ] && awk 'NR==1{print $1}' $MC 2>/dev/null || true)
  [ -n "$COST" ] && echo "    bbregister min cost $COST (FreeSurfer: < ~0.6 is typical for a good fit)"
else
  echo "[b] SynthSeg --parc --robust on the T1"
  mri_synthseg --i $B/data/${S}_t1.nii.gz --o $W/${S}_t1synthseg.nii.gz --parc --robust --threads ${CSS_THREADS:-8} > /dev/null || die "mri_synthseg failed"
  LAB=$W/${S}_t1synthseg.nii.gz
  echo "[c] register SWI -> T1 (mri_coreg, rigid, mutual information)"
  mri_coreg --mov $SWI --ref $W/${S}_t1_brain.nii.gz --reg $W/${S}_swi2t1.lta --dof 6 > $W/${S}_coreg.log 2>&1 \
    || die "mri_coreg failed (see $W/${S}_coreg.log)"
  REGM=mri_coreg; COST=""
fi

# v4.14 QC: registration method + cost (bbregister mincost; mri_coreg reports none) -> review/ID_qc.json
echo "$COST" | grep -Eq '^[0-9.eE+-]+$' || COST=null
printf '{"method": "%s", "cost": %s}\n' "$REGM" "$COST" > $W/${S}_reg.json
echo "[d] labels -> SWI grid (nearest neighbour)"
mri_vol2vol --mov $SWI --targ $LAB --lta $W/${S}_swi2t1.lta --inv --interp nearest \
  --o $W/${S}_t1seg_swispace.nii.gz > /dev/null || die "mri_vol2vol (labels -> SWI) failed"
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
echo "QC - labels must hug the SWI cortex in the CORONAL and SAGITTAL views too (an up-down slip is"
echo "invisible in axial slices); FLAIR must line up (toggle layers). detect_css.py reports label_shift_mm:"
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
    grep -m2 -E "raw components|anatomy " $B/logs/${s}_run.log
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
# v4.18.3: copy to a temporary name, then rename - a batch that is running keeps reading its old file
# (cp over a running zsh script changes it under the shell's feet)
for f in $SRC/*.py $SRC/*.sh; do
  t=$B/scripts/.$(basename $f).new; cp $f $t && chmod +x $t && mv -f $t $B/scripts/$(basename $f) || { echo "install FAILED at $f"; exit 1; }
done
grep -q 'css_project/scripts' ~/.zshrc || echo 'export PATH="$HOME/css_project/scripts:$PATH"' >> ~/.zshrc
echo "installed to $B/scripts  - open a new Terminal window"
```

## scripts/agreement.py

```python
#!/usr/bin/env python3
"""Reproducibility statistics as in van Harten et al. (NeuroImage Clin 2023;38:103447).
usage: agreement.py --a DIR_A --b DIR_B [--col grown_volume_mm3] [--plot out.png]
       agreement.py --vs-score DIR [--col surface_mm2] [--plot out.png]
  DIR = a copy of the review folder after one complete scoring session, e.g.
        cp -R ~/css_project/review ~/css_project/review_reader1_session1
        (it must contain css_scores.csv and, for Dice, ID_css_mask.nii.gz)
  --a/--b      two sessions of one reader (intra-observer) or two readers (inter-observer):
               Pearson r, ICC(A,1) absolute agreement (two-way), Bland-Altman bias and 95 % limits of
               agreement, and per-subject Dice of the cSS masks. van Harten: ICC 0.995, Pearson 0.991,
               mean Dice 0.75 (Dice is low for small, sparse masks even when volumes agree).
               v4.17: quadratic weighted kappa for multifocality_0_4 and L/R_score (ordinal visual scale).
  --vs-score   size measure against the 0-4 multifocality score (their Fig. 4: wide spread inside the
               top categories = the ceiling effect a continuous measure avoids)
Compare sizes only between scans acquired with the same protocol (field strength, TE, voxel size)."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, numpy as np, pandas as pd, nibabel as nib

arg = lambda f, d="": sys.argv[sys.argv.index(f) + 1] if f in sys.argv else d
col = arg("--col", "grown_volume_mm3")
ex = lambda p: os.path.expanduser(p)


def icc_a1(x, y):
    """two-way, absolute agreement, single measures (McGraw & Wong ICC(A,1))"""
    Y = np.c_[x, y]; n, k = Y.shape
    gm = Y.mean(); rm = Y.mean(1); cm = Y.mean(0)
    msr = k * ((rm - gm) ** 2).sum() / (n - 1)
    msc = n * ((cm - gm) ** 2).sum() / (k - 1)
    mse = ((Y - rm[:, None] - cm[None, :] + gm) ** 2).sum() / ((n - 1) * (k - 1))
    return float((msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n))


if arg("--vs-score"):
    d = pd.read_csv(f"{ex(arg('--vs-score'))}/css_scores.csv")
    if col not in d.columns: sys.exit(f"column {col} not in css_scores.csv (re-score with score_css.py v4.13)")
    print(f"{col} by multifocality score ({len(d)} subjects):")
    for s, g in d.groupby("multifocality_0_4"):
        print(f"  score {s}: n={len(g)}  median {g[col].median():.1f}  range {g[col].min():.1f}-{g[col].max():.1f}")
    if arg("--plot"):
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 4))
        jit = np.random.default_rng(0).uniform(-0.12, 0.12, len(d))
        ax.scatter(d.multifocality_0_4 + jit, d[col], s=25)
        ax.set_xlabel("cSS multifocality score (0-4)"); ax.set_ylabel(col); ax.set_xticks(range(5))
        fig.tight_layout(); fig.savefig(ex(arg("--plot")), dpi=150); print(f"plot -> {arg('--plot')}")
    sys.exit(0)

if not (arg("--a") and arg("--b")): sys.exit(__doc__)
A, Bd = ex(arg("--a")), ex(arg("--b"))
da, db = pd.read_csv(f"{A}/css_scores.csv"), pd.read_csv(f"{Bd}/css_scores.csv")
m = da[["subject", col]].merge(db[["subject", col]], on="subject", suffixes=("_a", "_b")).dropna()
if len(m) < 2: sys.exit(f"need >=2 subjects scored in both folders (found {len(m)})")
x, y = m[f"{col}_a"].to_numpy(float), m[f"{col}_b"].to_numpy(float)
diff = x - y; bias = diff.mean(); sd = diff.std(ddof=1)
r = float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else float("nan")
print(f"{len(m)} subjects, {col}:  Pearson r {r:.3f}   ICC(A,1) {icc_a1(x, y):.3f}   "
      f"Bland-Altman bias {bias:.1f}, 95 % limits {bias - 1.96 * sd:.1f} to {bias + 1.96 * sd:.1f}")
# v4.17: weighted kappa for the ordinal multifocality score (and per hemisphere), as used for visual rating scales
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import cohen_kappa
for sc, cats in (("multifocality_0_4", range(5)), ("L_score", range(3)), ("R_score", range(3))):
    if sc in da.columns and sc in db.columns:
        k_ = da[["subject", sc]].merge(db[["subject", sc]], on="subject", suffixes=("_a", "_b")).dropna()
        if len(k_) >= 2:
            print(f"{sc}: quadratic weighted kappa {cohen_kappa(k_[sc + '_a'].astype(int), k_[sc + '_b'].astype(int), cats, 'quadratic'):.3f}"
                  f"  exact agreement {(k_[sc + '_a'] == k_[sc + '_b']).mean():.0%}  (n={len(k_)})")
dice = []
for s in m.subject:
    pa, pb = f"{A}/{s}_css_mask.nii.gz", f"{Bd}/{s}_css_mask.nii.gz"
    if os.path.exists(pa) and os.path.exists(pb):
        a_, b_ = nib.load(pa).get_fdata() > 0, nib.load(pb).get_fdata() > 0
        if a_.shape == b_.shape and (a_.any() or b_.any()):
            dice.append((s, 2 * (a_ & b_).sum() / (a_.sum() + b_.sum())))
if dice:
    dv = np.array([v for _, v in dice])
    print(f"Dice of the cSS masks: mean {dv.mean():.2f} +/- {dv.std():.2f}  ("
          + ", ".join(f"{s} {v:.2f}" for s, v in dice) + ")")
print(m.to_string(index=False))
if arg("--plot"):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    lim = [0, max(x.max(), y.max()) * 1.05]
    ax[0].scatter(x, y); ax[0].plot(lim, lim, "k--", lw=0.8)
    ax[0].set_xlabel(f"{col} (A)"); ax[0].set_ylabel(f"{col} (B)"); ax[0].set_title(f"r = {r:.3f}")
    mean = (x + y) / 2; ax[1].scatter(mean, diff)
    for v, ls in ((bias, "-"), (bias - 1.96 * sd, "--"), (bias + 1.96 * sd, "--")): ax[1].axhline(v, color="k", ls=ls, lw=0.8)
    ax[1].set_xlabel("mean of A and B"); ax[1].set_ylabel("A - B"); ax[1].set_title("Bland-Altman")
    fig.tight_layout(); fig.savefig(ex(arg("--plot")), dpi=150); print(f"plot -> {arg('--plot')}")
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

## scripts/check_known.py

```python
#!/usr/bin/env python3
"""Safety check: are KNOWN lesions still detected after a detector change?
usage: check_known.py SUBJECT --old OLD_candidates.nii.gz --ids 3,4,5
For each listed lesion of an older candidate map (e.g. review/v3/P006_candidates.nii.gz) it reports
the current candidate(s) covering it, or the exclusion reason (review/ID_excluded.csv), or
"NOT DETECTED". Grids must match (same imported SWI)."""
import sys, os, numpy as np, pandas as pd, nibabel as nib
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base

subj = sys.argv[1]; base = _base()
arg = lambda f: sys.argv[sys.argv.index(f) + 1] if f in sys.argv else ""
if not arg("--old") or not arg("--ids"): sys.exit(__doc__)
old = nib.load(os.path.expanduser(arg("--old"))).get_fdata().astype(int)
cand = nib.load(f"{base}/review/{subj}_candidates.nii.gz").get_fdata().astype(int)
xp = f"{base}/review/{subj}_excluded.nii.gz"
excl = nib.load(xp).get_fdata().astype(int) if os.path.exists(xp) else np.zeros_like(cand)
xr = pd.read_csv(f"{base}/review/{subj}_excluded.csv").set_index("excl_id")["reason"] \
    if os.path.exists(f"{base}/review/{subj}_excluded.csv") else pd.Series(dtype=str)
cd = pd.read_csv(f"{base}/review/{subj}_candidates.csv").set_index("cand_id")
def info(v):     # v4.18 recorded-only extras: review tier and what the ends do
    if v not in cd.index or "review_tier" not in cd.columns: return f"#{v}"
    r = cd.loc[v]; t = f"#{v} (tier {int(r.review_tier)}"
    t += f": {r.tier_reason}" if isinstance(r.tier_reason, str) and r.tier_reason else ""
    return t + f"; ends {r.end_type_max}/{r.end_type_min})" if isinstance(r.end_type_max, str) else t + ")"
if old.shape != cand.shape: sys.exit("grids differ - re-import the same SWI")
kept = 0
for i in [int(x) for x in arg("--ids").split(",") if x.strip()]:
    m = ndi.binary_dilation(old == i, iterations=2)
    c = sorted(set(np.unique(cand[m])) - {0}); x = sorted(set(np.unique(excl[m])) - {0})
    if c: kept += 1
    msg = (f"kept as candidate {', '.join(info(v) for v in c)}" if c else "") + \
          ("; " if c and x else "") + ("; ".join(f"EXCLUDED ({xr.get(v, '?')})" for v in x) if x else "")
    print(f"known lesion #{i}: {msg or 'NOT DETECTED'}")
print(f"{subj}: {kept} of {len(arg('--ids').split(','))} known lesions still shown for review")
```

## scripts/consensus.py

```python
#!/usr/bin/env python3
"""Two-reader consensus, adjudication and agreement statistics for expert candidate labels (v4.17).
usage: consensus.py --readers expertA,expertB [--subjects P012,P013] [--session 1]
       consensus.py --readers expertA,expertB --apply labels/adjudication_todo.csv
Reads labels/candidate_labels.csv (expert_import.py). For every object both readers labelled:
  same fine label            -> consensus = that label                     (method "agreement")
  same class, different code -> consensus = the class (e.g. other_vein)    (method "agreement")
  different class            -> listed in labels/adjudication_todo.csv     (method "pending")
Classes (css_common.label_class): cSS, uncertain, surface_vein, other_vein, artifact, other_normal.
Objects only one reader saw are kept as "single_reader". A third reader (or a consensus meeting) fills
adjudicated_label (a sheet code such as C / VS / A, or a label name) and adjudicator, then --apply writes them
as method "adjudicated". Output: labels/consensus.csv (one row per object) and labels/agreement_stats.json.
Agreement: Cohen's kappa for cSS vs not-cSS and for the six classes, raw agreement, with 95 % CIs from
bootstrapping PATIENTS (candidates of one patient are not independent)."""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, LABELS, label_class, cohen_kappa, append_rows

arg = lambda f, d="": sys.argv[sys.argv.index(f) + 1] if f in sys.argv else d
if not arg("--readers") or len(arg("--readers").split(",")) != 2: sys.exit(__doc__)
ra, rb = arg("--readers").split(","); base = _base()
lp = f"{base}/labels/candidate_labels.csv"
if not os.path.exists(lp): sys.exit(f"no {lp} - import expert answers first (expert_import.py)")
L = pd.read_csv(lp)
session = int(arg("--session", "1")); L = L[L.session == session]
if arg("--subjects"): L = L[L.study_id.astype(str).isin(arg("--subjects").split(","))]
CLASSES = ["cSS", "uncertain", "surface_vein", "other_vein", "artifact", "other_normal"]
cons_p = f"{base}/labels/consensus.csv"; todo_p = f"{base}/labels/adjudication_todo.csv"
key = ["study_id", "cand_uid", "list"]

if arg("--apply"):
    A = pd.read_csv(os.path.expanduser(arg("--apply")))
    A = A[A.adjudicated_label.notna() & (A.adjudicated_label.astype(str).str.strip() != "")]
    names = {v[0] for v in LABELS.values()}
    rows = []
    for _, r in A.iterrows():
        v = str(r.adjudicated_label).strip()
        lab = LABELS[v.upper()][0] if v.upper() in LABELS else v if v in names else None
        if lab is None: sys.exit(f"cannot read adjudicated_label '{v}' (use a sheet code or a label name)")
        rows.append(dict(study_id=r.study_id, cand_uid=r.cand_uid, list=r.list, obj_id=r.obj_id,
                         consensus_label=lab, consensus_class=label_class(lab), method="adjudicated",
                         readers=f"{ra}+{rb}", adjudicator=r.get("adjudicator", "")))
    append_rows(cons_p, rows, key)
    print(f"{len(rows)} adjudicated labels -> {cons_p}"); sys.exit(0)

a = L[L.reader_id == ra].set_index(key); b = L[L.reader_id == rb].set_index(key)
if not len(a) or not len(b): sys.exit(f"need labels from both {ra} and {rb} (session {session})")
both = a.index.intersection(b.index)
rows, todo = [], []
for k in a.index.union(b.index):
    if k in both:
        la, lb = a.loc[k, "label"], b.loc[k, "label"]; ca, cb = label_class(la), label_class(lb)
        r = dict(zip(key, k), obj_id=a.loc[k, "obj_id"], readers=f"{ra}+{rb}", adjudicator="")
        if ca == cb:
            rows.append(dict(r, consensus_label=la if la == lb else ca, consensus_class=ca, method="agreement"))
        else:
            rows.append(dict(r, consensus_label="", consensus_class="", method="pending"))
            todo.append(dict(zip(key, k), obj_id=a.loc[k, "obj_id"], sheet_no=a.loc[k, "sheet_no"],
                             **{f"label_{ra}": la, f"label_{rb}": lb}, adjudicated_label="", adjudicator=""))
    else:
        src = a if k in a.index else b
        rows.append(dict(zip(key, k), obj_id=src.loc[k, "obj_id"], consensus_label=src.loc[k, "label"],
                         consensus_class=label_class(src.loc[k, "label"]), method="single_reader",
                         readers=ra if k in a.index else rb, adjudicator=""))
append_rows(cons_p, rows, key)
pd.DataFrame(todo, columns=key + ["obj_id", "sheet_no", f"label_{ra}", f"label_{rb}", "adjudicated_label",
                                  "adjudicator"]).to_csv(todo_p, index=False)

# agreement on objects both readers labelled; CIs by bootstrapping patients
P = pd.DataFrame([dict(study=k[0], ca=label_class(a.loc[k, "label"]), cb=label_class(b.loc[k, "label"])) for k in both])
def stats(D):
    return dict(kappa_css=cohen_kappa(D.ca == "cSS", D.cb == "cSS", [False, True]),
                kappa_class=cohen_kappa(D.ca, D.cb, CLASSES), agreement=float((D.ca == D.cb).mean()))
est = stats(P); rng = np.random.default_rng(0); studies = P.study.unique(); boot = []
if len(studies) >= 3:
    for _ in range(1000):
        pick = rng.choice(studies, len(studies), replace=True)
        boot.append(stats(pd.concat([P[P.study == s_] for s_ in pick])))
ci = {k: ([round(float(np.nanpercentile([b_[k] for b_ in boot], q)), 3) for q in (2.5, 97.5)] if boot else None)
      for k in est}
out = dict(readers=[ra, rb], session=session, n_objects=int(len(P)), n_patients=int(len(studies)),
           estimates={k: round(v, 3) for k, v in est.items()}, ci95_patient_bootstrap=ci,
           confusion=pd.crosstab(P.ca, P.cb).to_dict(), pending=len(todo))
json.dump(out, open(f"{base}/labels/agreement_stats.json", "w"), indent=2)
f = lambda k: f"{est[k]:.2f}" + (f" [{ci[k][0]:.2f}, {ci[k][1]:.2f}]" if ci[k] else " (CI needs >= 3 patients)")
print(f"{ra} vs {rb}: {len(P)} objects in {len(studies)} patient(s): kappa cSS-vs-not {f('kappa_css')}, "
      f"kappa 6 classes {f('kappa_class')}, raw agreement {f('agreement')}")
print(f"consensus -> {cons_p} ({sum(r['method'] == 'agreement' for r in rows)} agreed, {len(todo)} pending, "
      f"{sum(r['method'] == 'single_reader' for r in rows)} single-reader)"
      + (f"; to adjudicate: {todo_p} (fill adjudicated_label, then --apply)" if todo else ""))
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


# ---- v4.15 expert label vocabulary (docs/AUDIT_v4.13_and_FP_plan.md section 4)
# sheet code -> (fine label, coarse call). The coarse call is what review_css.py / feature_report.py /
# mark_css.py already understand, so old single letters (C V N H U) keep their meaning.
LABELS = {
    "C":   ("cSS", "cSS"),
    "V":   ("vein_unspecified", "Vein"),
    "VS":  ("vein_surface", "Vein"),            # cortical vein running ALONG the pial surface (hard negative)
    "VC":  ("vein_sulcal", "Vein"),             # vein in the middle of the sulcal CSF
    "VT":  ("vein_transcortical", "Vein"),      # medullary / transcortical, perpendicular to the cortex
    "VX":  ("vein_sinus_adjacent", "Vein"),     # vein joining / next to a dural sinus
    "TV":  ("thrombosed_vein", "Vein"),
    "A":   ("artifact_airbone", "Artifact"),    # skull base / orbitofrontal / temporal pole susceptibility
    "AM":  ("artifact_motion", "Artifact"),     # motion, ghosting, Gibbs, registration
    "N":   ("normal_other", "Normal"),
    "NI":  ("normal_iron_cortex", "Normal"),    # e.g. iron-rich motor cortex
    "MB":  ("microbleed", "Normal"),
    "CA":  ("calcification", "Normal"),
    "SAH": ("acute_cSAH", "Normal"),
    "LN":  ("laminar_necrosis", "Normal"),
    "HI":  ("haemorrhagic_infarct", "Normal"),
    "IS":  ("infratentorial_SS", "Normal"),
    "H":   ("ICH_related", "Near ICH"),
    "U":   ("uncertain", "Unsure"),
}
LABEL_HELP = [("C", "cortical superficial siderosis"),
              ("V / VS / VC / VT / VX", "vein: any / along surface / mid-sulcus / transcortical / at sinus"),
              ("TV", "thrombosed vein"),
              ("A / AM", "artifact: air-bone susceptibility / motion, ghosting"),
              ("N / NI", "normal or other / normal iron-rich cortex"),
              ("MB / CA", "microbleed / calcification"),
              ("SAH", "acute convexity SAH (FLAIR-bright sulcus)"),
              ("LN / HI", "laminar necrosis / haemorrhagic infarct"),
              ("IS", "infratentorial (classical) superficial siderosis"),
              ("H", "siderosis contiguous with a lobar ICH"),
              ("U", "uncertain")]


def parse_answer(tok):
    """'12VS' or '12VS:4' -> (12, code, confidence or None); raises ValueError"""
    import re
    m = re.fullmatch(r"(\d+)([A-Z]+)(?::([1-5]))?", tok.strip().upper())
    if not m or m.group(2) not in LABELS:
        raise ValueError(f"cannot read '{tok}' (e.g. 12C, 12VS, 12VS:4; codes {', '.join(LABELS)})")
    return int(m.group(1)), m.group(2), (int(m.group(3)) if m.group(3) else None)


def cand_uid(subj, r):
    """stable candidate identity across detector re-runs (same fingerprint as match_calls)"""
    return f"{subj}:{int(r['vox_i'])}:{int(r['vox_j'])}:{int(r['vox_k'])}:{float(r['volume_mm3']):g}"


def append_rows(path, rows, key):
    """append rows to a CSV table, replacing earlier rows with the same key columns"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    new = pd.DataFrame(rows)
    if os.path.exists(path):
        old = pd.read_csv(path)
        if len(new):
            k_old = old[key].astype(str).agg("|".join, axis=1)
            k_new = set(new[key].astype(str).agg("|".join, axis=1))
            old = old[~k_old.isin(k_new)]
        new = pd.concat([old, new], ignore_index=True)
    new.to_csv(path, index=False)
    return new


# ---- v4.17 provenance helpers
def fresh(path, ref):
    """True if path exists and is not older than ref (a derived file made after its source was imported)"""
    return os.path.exists(path) and (not os.path.exists(ref) or os.path.getmtime(path) >= os.path.getmtime(ref))


def sha256(path, block=1 << 20):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def update_manifest(b, subj, section, data):
    """merge {section: data} into review/ID_manifest.json (one machine-readable provenance record per case)"""
    import json
    p = f"{b}/review/{subj}_manifest.json"
    m = json.load(open(p)) if os.path.exists(p) else {"subject": subj}
    m[section] = data
    os.makedirs(os.path.dirname(p), exist_ok=True)
    json.dump(m, open(p, "w"), indent=2, default=str)


def environment_record():
    """python / package / FreeSurfer versions for the manifest"""
    import platform, importlib
    rec = {"python": platform.python_version(), "platform": platform.platform()}
    for mod in ("numpy", "scipy", "nibabel", "pandas", "skimage", "matplotlib"):
        try:
            rec[mod] = importlib.import_module(mod).__version__
        except Exception:
            rec[mod] = None
    fs = os.environ.get("FREESURFER_HOME", "")
    stamp = os.path.join(fs, "build-stamp.txt")
    rec["freesurfer"] = open(stamp).read().strip() if os.path.exists(stamp) else (fs or None)
    return rec


def label_class(label):
    """fine expert label -> display / agreement class (shared by surface_overlay.py and consensus.py)"""
    if label == "cSS": return "cSS"
    if label == "uncertain": return "uncertain"
    if label == "vein_surface": return "surface_vein"
    if label.startswith("vein_") or label == "thrombosed_vein": return "other_vein"
    if label.startswith("artifact"): return "artifact"
    return "other_normal"


def cohen_kappa(a, b, categories=None, weights=None):
    """Cohen's kappa for two raters; weights=None (unweighted), 'linear' or 'quadratic' (ordered categories)"""
    a, b = list(a), list(b)
    cats = list(categories) if categories is not None else sorted(set(a) | set(b))
    if len(a) == 0 or len(cats) < 2:
        return float("nan")
    k = {c: i for i, c in enumerate(cats)}; n = len(cats)
    O = np.zeros((n, n))
    for x, y in zip(a, b):
        O[k[x], k[y]] += 1
    O /= O.sum()
    E = np.outer(O.sum(1), O.sum(0))
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    W = np.where(i == j, 0.0, 1.0) if weights is None else (np.abs(i - j) / (n - 1)) ** (2 if weights == "quadratic" else 1)
    den = (W * E).sum()
    return float(1 - (W * O).sum() / den) if den > 0 else float("nan")


def skeleton_length(sk, vox):
    """v4.18: centre-line length (mm) of a 2-D skeleton from its steps (detect_css used pixels x vox[0], which
    undercounts diagonals): orthogonal steps vox[0] / vox[1], diagonal steps their hypotenuse; a diagonal is
    skipped where an orthogonal corner already joins the two pixels"""
    import numpy as np
    sk = np.asarray(sk, bool); a, b = float(vox[0]), float(vox[1]); d = float(np.hypot(a, b))
    h = sk[1:, :] & sk[:-1, :]; v = sk[:, 1:] & sk[:, :-1]
    d1 = sk[1:, 1:] & sk[:-1, :-1] & ~(sk[1:, :-1] | sk[:-1, 1:])
    d2 = sk[1:, :-1] & sk[:-1, 1:] & ~(sk[1:, 1:] | sk[:-1, :-1])
    return float(h.sum() * a + v.sum() * b + (d1.sum() + d2.sum()) * d)


def radius_profile(mask, vox, sample=None):
    """v4.18: in-plane half-width (mm) along the 2-D centre line, slice by slice (cf. Avadiappan 2020: radius from a
    Euclidean distance transform). Distance to the nearest outside voxel (sampling vox[0], vox[1]) minus half a voxel,
    read on the skeleton of `sample` (default: the mask itself), so a 1-voxel-wide line reads half a voxel and a
    3-voxel-wide band 1.5 voxels. Quantised at the in-plane voxel size: compare only within one protocol."""
    import numpy as np
    from scipy import ndimage as ndi
    from skimage.morphology import skeletonize
    mask = np.asarray(mask, bool); sample = mask if sample is None else np.asarray(sample, bool)
    half = 0.5 * float(np.sqrt(vox[0] * vox[1])); out = []
    for k in range(mask.shape[2]):
        s = sample[:, :, k] & mask[:, :, k]
        if not s.any(): continue
        e = ndi.distance_transform_edt(np.pad(mask[:, :, k], 1), sampling=vox[:2])[1:-1, 1:-1]
        sk = skeletonize(s) & s
        out.append(e[sk] - half)
    return np.concatenate(out) if out else np.zeros(0)


GAMMA_HZ_T = 42.577e6       # proton gyromagnetic ratio / 2 pi


def susceptibility_signal(chi_ppm, vox, affine, b0_t=3.0, te_s=0.020, pad=16):
    """v4.18 susceptibility forward model for synthetic lesions (cf. CenSynCMB 2026: dipole kernel + static dephasing).
    chi_ppm: susceptibility DIFFERENCE map (ppm, 0 = host tissue). The field perturbation (ppm of B0) is the
    convolution with the unit dipole, done in k-space: D(k) = 1/3 - (k.b0)^2 / |k|^2, with b0 the scanner z axis
    expressed in the array axes (oblique acquisitions handled through the affine). Phase = 2 pi gamma B0 TE field.
    Intravoxel dephasing for a field that varies linearly across the voxel: |sinc| along each axis of the local
    phase gradient x voxel size. Returns (phase_rad, dephasing_attenuation 0-1). An approximation for
    augmentation and stress tests - not a full GRE / SWI / QSM simulator."""
    import numpy as np
    chi = np.pad(np.asarray(chi_ppm, np.float64), pad)
    vox = np.asarray(vox, float)
    R = np.asarray(affine, float)[:3, :3] / vox[None, :]
    b0 = R.T @ np.array([0.0, 0.0, 1.0]); b0 /= np.linalg.norm(b0)          # B0 (world z) in array axes
    k = np.meshgrid(*[np.fft.fftfreq(n, d=v) for n, v in zip(chi.shape, vox)], indexing="ij")
    k2 = k[0] ** 2 + k[1] ** 2 + k[2] ** 2
    kb = k[0] * b0[0] + k[1] * b0[1] + k[2] * b0[2]
    with np.errstate(invalid="ignore", divide="ignore"):
        D = 1.0 / 3.0 - np.where(k2 > 0, kb ** 2 / k2, 0.0)
    D[0, 0, 0] = 0.0
    field = np.real(np.fft.ifftn(np.fft.fftn(chi) * D))[pad:-pad, pad:-pad, pad:-pad]       # ppm
    phase = 2 * np.pi * GAMMA_HZ_T * b0_t * te_s * field * 1e-6
    att = np.ones_like(phase)
    for ax in range(3):
        g = np.gradient(phase, vox[ax], axis=ax)                                # rad / mm
        att *= np.abs(np.sinc(g * vox[ax] / (2 * np.pi)))
    return phase.astype(np.float32), att.astype(np.float32)
```

## scripts/css_qc.py

```python
"""Per-case quality control for detect_css.py (v4.14) -> review/ID_qc.json.

Called by detect_css.py with the arrays it already has. Records the acquisition, orientation, coverage,
anatomy source and how well the label (pial) boundary agrees with the cortex/CSF edge seen on the SWI,
and sets a case status:
  ok                     nothing below triggered
  low_confidence_review  anatomy-based features are less trustworthy (SWI-only anatomy, poor edge
                         agreement, high registration cost, thick slices)
  insufficient_quality   the detector's assumptions do not hold (array axis 2 is not the slice /
                         superior-inferior axis, or labels cover too little of the brain)
v4.14: REPORTED ONLY. The status and the thresholds below are provisional (not validated on expert-
labelled cases) and do not change any candidate, exclusion or score."""
import os, json, glob, numpy as np, nibabel as nib
from scipy import ndimage as ndi

EDGE_LOW_CONF_MM = 1.5     # median |label pial boundary - SWI cortex/CSF edge| above this -> low confidence
REG_COST_LOW_CONF = 0.8    # bbregister mincost above this -> low confidence (FreeSurfer typical < ~0.6)
LABEL_SLICE_MIN = 0.9      # fraction of brain-containing slices that also contain labels
THICK_SLICE_MM = 2.0       # slices thicker than this -> low confidence (sub-voxel pial thresholds)
EDGE_WIN_MM = 3.0          # search for the SWI edge within +-3 mm of the label boundary, along the normal


def read_acq(base, subj):
    """acquisition parameters from the dcm2niix JSON of the (non-phase) SWI series, if present"""
    for jp in sorted(glob.glob(f"{base}/raw/{subj}/SWI/*.json")):
        try:
            js = json.load(open(jp))
        except Exception:
            continue
        it = " ".join(js.get("ImageType", [])).upper()
        if "PHASE" in it or " P " in f" {it} ":
            continue
        te = js.get("EchoTime")
        return dict(sequence="SWI" if "SWI" in it else ("T2*-GRE" if it else "unknown"),
                    image_type=it, manufacturer=js.get("Manufacturer", ""),
                    field_T=js.get("MagneticFieldStrength"),
                    TE_ms=round(te * 1000, 1) if isinstance(te, (int, float)) else None,
                    slice_thickness_mm=js.get("SliceThickness"), acq_dim=js.get("MRAcquisitionType"))
    return {}


def edge_agreement(If, brain, pial_sd, normal, vox, sel, rng, normals=False):
    """Signed offset (mm) of the strongest dark->bright step (cortex -> CSF) along the outward normal,
    searched within +-EDGE_WIN_MM of label-boundary voxels. 0 = the labels' pial boundary sits on the
    SWI edge. Returns offsets of up to 4000 sampled boundary voxels in `sel` (and, with normals=True, the
    outward unit normals of those voxels, for the shift fit)."""
    pts = np.argwhere(sel)
    if len(pts) < 50:
        return (np.array([]), np.zeros((0, 3))) if normals else np.array([])
    if len(pts) > 4000:
        pts = pts[rng.choice(len(pts), 4000, replace=False)]
    n = np.stack([normal[a][tuple(pts.T)] for a in range(3)], 1)
    t = np.arange(-EDGE_WIN_MM, EDGE_WIN_MM + 1e-6, 0.25, dtype=np.float32)
    q = pts[:, None, :] + t[None, :, None] * n[:, None, :] / np.array(vox)
    co = q.reshape(-1, 3).T
    prof = ndi.map_coordinates(If, co, order=1, mode="nearest").reshape(q.shape[:2])
    inb = ndi.map_coordinates(brain.astype(np.float32), co, order=0, mode="constant").reshape(q.shape[:2]) > 0
    d = np.diff(prof, axis=1)
    d[~(inb[:, 1:] & inb[:, :-1])] = -np.inf          # steps that leave the brain mask do not count
    j = np.argmax(d, axis=1); ok = np.isfinite(d[np.arange(len(d)), j]) & (d[np.arange(len(d)), j] > 0)
    off = (0.5 * (t[j] + t[j + 1]))[ok]
    return (off, n[ok]) if normals else off


def _offsets(If, brain, normal, vox, pts, centre):
    """offset (mm, along the outward normal) of the strongest dark->bright step within centre +- EDGE_WIN_MM"""
    n = np.stack([normal[a][tuple(pts.T)] for a in range(3)], 1)
    t = np.arange(-EDGE_WIN_MM, EDGE_WIN_MM + 1e-6, 0.25, dtype=np.float32)
    tt = centre[:, None] + t[None, :]
    q = pts[:, None, :] + tt[:, :, None] * n[:, None, :] / np.array(vox)
    co = q.reshape(-1, 3).T
    prof = ndi.map_coordinates(If, co, order=1, mode="nearest").reshape(q.shape[:2])
    inb = ndi.map_coordinates(brain.astype(np.float32), co, order=0, mode="constant").reshape(q.shape[:2]) > 0
    d = np.diff(prof, axis=1)
    d[~(inb[:, 1:] & inb[:, :-1])] = -np.inf
    j = np.argmax(d, axis=1); r = np.arange(len(d))
    ok = np.isfinite(d[r, j]) & (d[r, j] > 0)
    off = np.where(ok, 0.5 * (tt[r, j] + tt[r, j + 1]), np.nan)
    return off, n


def shift_fit(If, brain, normal, vox, sel, rng, iters=6):
    """Rigid label->SWI shift s (mm): a label map displaced by s moves its boundary by s, so the SWI edge appears
    at -s.n along each outward normal n. Iterative: each search window is re-centred on the current estimate
    (a fixed +-3 mm window cannot see a large shift), offsets > 1.5 mm from the fit are dropped as outliers
    (veins, the white-matter edge). Returns (s, bias, scatter, share of samples with an edge): bias = median
    offset left after the shift (labels too small > 0 / too large < 0), scatter = median absolute residual."""
    pts = np.argwhere(sel)
    if len(pts) < 50:
        return None, None, None, None
    if len(pts) > 4000:
        pts = pts[rng.choice(len(pts), 4000, replace=False)]
    s = np.zeros(3)
    for k in range(iters):
        n = np.stack([normal[a][tuple(pts.T)] for a in range(3)], 1)
        off, n = _offsets(If, brain, normal, vox, pts, -(n @ s))
        ok = np.isfinite(off)
        if ok.sum() < 50:
            return None, None, None, float(ok.mean())
        res = off + n @ s
        if k:                                          # robust: drop outliers around the current fit
            ok &= np.abs(res - np.median(res[ok])) < 1.5
        s_new, *_ = np.linalg.lstsq(n[ok], -off[ok], rcond=None)
        done = np.linalg.norm(s_new - s) < 0.05
        s = s_new
        if done:
            break
    res = (off + n @ s)[np.isfinite(off)]
    return s, float(np.median(res)), float(np.median(np.abs(res - np.median(res)))), float(np.isfinite(off).mean())


def case_qc(base, subj, swi, I, If, seg, brain, tissue, cortex, pial_sd, normal, rim, hemiL, vox,
            t1_anat, full_bottom):
    rng = np.random.default_rng(0)
    qc = dict(subject=subj, qc_version="v4.19.3")
    qc["acquisition"] = read_acq(base, subj)
    qc["voxel_mm"] = [round(v, 3) for v in vox]
    qc["orientation"] = "".join(nib.aff2axcodes(swi.affine))
    si_axis = int(np.argmax(np.abs(swi.affine[2, :3])))
    qc["si_axis"] = si_axis
    qc["thickest_axis"] = int(np.argmax(vox))
    # detect_css.py runs its 2-D ridge filter, slice counts and merge gap in array axis 2 and treats
    # axis 2 as superior-inferior ("long structure"); a sagittally/coronally stored SWI breaks that
    qc["orientation_ok"] = bool(si_axis == 2 and (vox[2] >= max(vox[0], vox[1]) - 1e-3))
    qc["coverage_includes_skull_base"] = bool(full_bottom)
    if os.path.exists(f"{base}/work/{subj}_a2009s_swispace.nii.gz") and t1_anat:
        qc["anatomy"] = "T1 recon-all"
    else:
        qc["anatomy"] = "T1" if t1_anat else "SWI-only"
    reg_p = f"{base}/work/{subj}_reg.json"
    qc["registration"] = json.load(open(reg_p)) if os.path.exists(reg_p) else None
    # label coverage: brain-containing slices (axis 2) that also contain labels; labelled share of brain
    bs = brain.any(axis=(0, 1)); ls = (seg > 0).any(axis=(0, 1))
    qc["label_slice_coverage"] = round(float((bs & ls).sum() / max(int(bs.sum()), 1)), 3)
    qc["label_voxel_frac"] = round(float(((seg > 0) & brain).sum() / max(int(brain.sum()), 1)), 3)
    # pial confidence: label boundary vs SWI cortex/CSF edge, overall and by hemisphere x height third
    # v4.18.3: the voxel layers on both sides of the label boundary sit at |pial_sd| >= the smallest voxel size,
    # so 0.75 x the largest voxel selected NOTHING on isotropic scans (LEMON 0.8 mm: edge_n 0, "not measurable")
    bnd = cortex & (np.abs(pial_sd) <= max(0.75 * max(vox), 1.01 * min(vox))) & ~rim
    off = edge_agreement(If, brain, pial_sd, normal, vox, bnd, rng)
    s_, bias_, spread_, found_ = shift_fit(If, brain, normal, vox, bnd, rng)
    # v4.17.1: separate a systematic misregistration (rigid shift, global bias) from scatter. On P006/P011 the
    # median |offset| (1.6 / 2.1 mm) mixed both; the status now uses the systematic parts only.
    qc["label_shift_mm"] = round(float(np.linalg.norm(s_)), 2) if s_ is not None else None
    qc["label_shift_vec_mm"] = [round(float(v), 2) for v in s_] if s_ is not None else None
    qc["edge_bias_mm"] = round(bias_, 2) if bias_ is not None else None
    qc["edge_scatter_mm"] = round(spread_, 2) if spread_ is not None else None
    qc["edge_found_frac"] = round(found_, 2) if found_ is not None else None
    qc["edge_offset_median_mm"] = round(float(np.median(off)), 2) if len(off) else None
    qc["edge_agreement_mm"] = round(float(np.median(np.abs(off))), 2) if len(off) else None
    qc["edge_n"] = int(len(off))
    zw = (swi.affine @ np.c_[np.argwhere(bnd), np.ones(int(bnd.sum()))].T)[2] if bnd.any() else np.array([])
    reg = {}
    if len(zw):
        lo, hi = np.percentile(zw, [0, 100]); cuts = np.linspace(lo, hi, 4)
        zmap = np.full(I.shape, -1, np.int8); zmap[bnd] = np.clip(np.searchsorted(cuts, zw, side="right") - 1, 0, 2)
        for h, hm in (("L", hemiL), ("R", ~hemiL)):
            for k, nm in enumerate(("inferior", "middle", "superior")):
                o = edge_agreement(If, brain, pial_sd, normal, vox, bnd & hm & (zmap == k), rng)
                reg[f"{h}_{nm}"] = round(float(np.median(np.abs(o))), 2) if len(o) >= 50 else None
    qc["edge_agreement_by_region_mm"] = reg
    # motion / ghosting proxy (informational): background noise of the un-stripped SWI vs brain signal
    raw_p = f"{base}/work/{subj}_raw.nii.gz"
    qc["background_ratio"] = None
    if os.path.exists(raw_p):
        raw = nib.load(raw_p)
        if raw.shape[:3] == I.shape:
            R = np.asarray(raw.dataobj, dtype=np.float32).reshape(I.shape)
            far = ndi.distance_transform_edt(~brain, sampling=vox) > 15.0
            if far.sum() > 1000 and brain.any():
                qc["background_ratio"] = round(float(np.std(R[far]) / (np.median(R[brain]) + 1e-6)), 4)
    # status (provisional thresholds, reported only)
    why_bad, why_low = [], []
    if not qc["orientation_ok"]:
        why_bad.append(f"array axis 2 is not the slice / S-I axis (orientation {qc['orientation']}, voxel {qc['voxel_mm']})")
    if qc["label_slice_coverage"] < LABEL_SLICE_MIN:
        why_bad.append(f"labels cover {qc['label_slice_coverage']:.0%} of brain slices")
    if qc["anatomy"] == "SWI-only":
        why_low.append("SWI-only anatomy (no T1)")
    if qc["edge_agreement_mm"] is None:
        why_low.append("pial edge agreement not measurable")
    else:
        if qc["label_shift_mm"] is not None and qc["label_shift_mm"] > EDGE_LOW_CONF_MM:
            why_low.append(f"labels shifted {qc['label_shift_mm']} mm against the SWI (rigid fit to the cortex edge)")
        if qc["edge_bias_mm"] is not None and abs(qc["edge_bias_mm"]) > EDGE_LOW_CONF_MM:
            why_low.append(f"label pial boundary {qc['edge_bias_mm']:+.2f} mm from the SWI cortex edge on average "
                           f"({'labels too small' if qc['edge_bias_mm'] > 0 else 'labels too large'})")
    cost = (qc["registration"] or {}).get("cost")
    # v4.19.3: the bbregister cost only lowers the status when the direct edge fit could not be measured. LEMON GRE
    # (P201/P202/P204): cost 0.95-0.96 every time while the labels sat 0.11-0.22 mm from the SWI cortex edge -
    # the T2*-weighted GRE contrast gives high costs even when aligned. The cost is still recorded.
    qc["notes"] = []
    if isinstance(cost, (int, float)) and cost > REG_COST_LOW_CONF:
        if qc["label_shift_mm"] is None:
            why_low.append(f"registration cost {cost}")
        else:
            qc["notes"].append(f"registration cost {cost} above {REG_COST_LOW_CONF} - not used, the edge fit was measured "
                               f"(label shift {qc['label_shift_mm']} mm)")
    if max(vox) > THICK_SLICE_MM:
        why_low.append(f"thick slices ({max(vox):.1f} mm)")
    qc["status"] = "insufficient_quality" if why_bad else ("low_confidence_review" if why_low else "ok")
    qc["reasons"] = why_bad + why_low
    qc["thresholds"] = dict(edge_low_conf_mm=EDGE_LOW_CONF_MM, reg_cost_low_conf=REG_COST_LOW_CONF,
                            label_slice_min=LABEL_SLICE_MIN, thick_slice_mm=THICK_SLICE_MM,
                            note="provisional, not validated; v4.14 reports only")
    os.makedirs(f"{base}/review", exist_ok=True)
    json.dump(qc, open(f"{base}/review/{subj}_qc.json", "w"), indent=2)
    return qc
```

## scripts/css_surface.py

```python
"""FreeSurfer pial surface -> SWI space, with a registration check (v4.16; shared by surface_features.py and
surface_overlay.py). FreeSurfer surfaces are in tkRAS. The transform is built from subjects/ID/mri/orig.mgz
(tkRAS -> T1 scanner RAS) and work/ID_swi2t1.lta (prep_anat.sh). The LTA's direction is NOT trusted: both
directions are tried and the one whose pial vertices lie on the label pial boundary of
work/ID_seg_swispace.nii.gz (which mri_vol2vol made with the same LTA) is used. surface_to_swi() stops the
script when even the better direction is more than AGREE_MM off, unless force=True."""
import os, sys, numpy as np, nibabel as nib
from scipy import ndimage as ndi
from nibabel.freesurfer.io import read_geometry

AGREE_MM = 2.0


def apply(T, X):
    return X @ T[:3, :3].T + T[:3, 3]


def read_lta(p):
    """4x4 RAS->RAS matrix of an LTA (type 1 RAS2RAS, or type 0 VOX2VOX converted with its volume info)"""
    L = open(p).read().splitlines()
    typ = int(next(l for l in L if l.strip().startswith("type")).split("=")[1].split()[0])
    i = next(k for k, l in enumerate(L) if l.strip().split()[:3] == ["1", "4", "4"])
    M = np.array([[float(x) for x in L[i + 1 + r].split()] for r in range(4)])
    if typ == 1:
        return M
    def v2r(tag):                                     # vox2ras from "src/dst volume info" block
        k = next(j for j, l in enumerate(L) if l.strip().startswith(tag))
        info = {}
        for l in L[k + 1:k + 10]:
            if "=" in l:
                key, val = l.split("=", 1); info[key.strip()] = val.split()
        dims = np.array([float(x) for x in info["volume"][:3]]); vs = np.array([float(x) for x in info["voxelsize"]])
        R = np.c_[[float(x) for x in info["xras"]], [float(x) for x in info["yras"]], [float(x) for x in info["zras"]]] * vs
        c = np.array([float(x) for x in info.get("cras", info.get("c_ras"))])   # FreeSurfer writes "cras"
        A = np.eye(4); A[:3, :3] = R; A[:3, 3] = c - R @ (dims / 2.0)
        return A
    return v2r("dst volume info") @ M @ np.linalg.inv(v2r("src volume info"))


def surface_to_swi(base, subj, sdir, surface_xfm="", force=False):
    """-> (T tkRAS->SWI world, source description, {option: median vertex-to-label-boundary mm}, surfs)
    surfs = {"lh": (vertices_tkRAS, faces), "rh": ...} of the pial surfaces"""
    swi = nib.load(f"{base}/data/{subj}_swi.nii"); vox = np.array(swi.header.get_zooms()[:3], float)
    seg = nib.load(f"{base}/work/{subj}_seg_swispace.nii.gz").get_fdata().astype(np.int32)
    brain = ndi.binary_fill_holes(np.asarray(swi.dataobj) > 0)
    tissue = ndi.binary_fill_holes(brain & (seg > 0) & (seg != 24))
    pial_sd = ndi.distance_transform_edt(~tissue, sampling=vox) - ndi.distance_transform_edt(tissue, sampling=vox)
    inv_aff = np.linalg.inv(swi.affine)
    agree = lambda Vw: float(np.median(np.abs(ndi.map_coordinates(pial_sd, apply(inv_aff, Vw).T, order=1, mode="nearest"))))
    surfs = {h: read_geometry(f"{sdir}/{subj}/surf/{h}.pial") for h in ("lh", "rh")}
    Vall = np.vstack([surfs["lh"][0], surfs["rh"][0]])[::50]
    if surface_xfm:
        T = np.loadtxt(surface_xfm).reshape(4, 4); src = "explicit --surface-xfm"
        checks = {"explicit": agree(apply(T, Vall))}
    else:
        orig = nib.load(f"{sdir}/{subj}/mri/orig.mgz")
        tkr2t1 = orig.header.get_vox2ras() @ np.linalg.inv(orig.header.get_vox2ras_tkr())
        lta = f"{base}/work/{subj}_swi2t1.lta"
        if os.path.getmtime(lta) < os.path.getmtime(f"{base}/data/{subj}_swi.nii") and not force:
            sys.exit(f"{subj}: {lta} is older than the SWI (re-imported?) - re-run prep_anat.sh (or --force)")
        R = read_lta(lta)
        opts = {"LTA maps SWI->T1": np.linalg.inv(R) @ tkr2t1, "LTA maps T1->SWI": R @ tkr2t1}
        checks = {k: agree(apply(Tk, Vall)) for k, Tk in opts.items()}
        src = min(checks, key=checks.get); T = opts[src]
    print(f"{subj}: surface -> SWI ({src}); median distance of pial vertices to the label pial boundary: "
          + ", ".join(f"{k} {v:.2f} mm" for k, v in checks.items()))
    if min(checks.values()) > AGREE_MM and not force:
        sys.exit(f"{subj}: surfaces do not match the SWI labels (best {min(checks.values()):.2f} mm > {AGREE_MM} mm)"
                 f" - check the registration / transform (or --force)")
    return T, src, checks, surfs
```

## scripts/css_vessel.py

```python
"""Bounded local vein tracking for detect_css.py (v4.16) - RECORDED ONLY, no rule uses these values.

vein_tree_mm saturated on real SWI because the dark tubular network merges brain-wide. These features
stay local and follow the vessel map (work/ID_vessels.nii.gz, the same VES map detect_css.py builds):

  cand_geodesic_mm      length of the candidate along its own voxels (two-sweep geodesic diameter)
  local_cont_max_mm     dark continuation beyond the candidate's ends, the longer of the two ends,
  local_cont_min_mm     ... and the shorter end. Each end is followed geodesically through vessel or dark
                        voxels (dark: so a vein draining into a wide dural sinus continues into it),
                        only forwards (inside a 60 deg cone around the end's outward tangent), up to
                        VES_REACH_MM; 1-voxel gaps are bridged at a cost
  local_cont_total_mm   sum of both ends
  terminal_ends         ends with < 3 mm continuation (2 = stops at both ends, cSS-like; 0 = runs on)
  local_net_mm          skeleton length of the vessel network connected to the candidate within reach
  local_branch_points   junctions in that network with >= 3 arms of >= 3 mm (one count per junction)
  local_branch_per10mm  junctions per 10 mm of network skeleton
  net_radius_med/p95_mm (v4.18) half-width of the connected vessel network (excluding the candidate), measured on
                        the dark + vessel voxels at the network's centre line (css_common.radius_profile)
  end_type_max/_min     (v4.18) what each end's continuation does, for the end with the longer / shorter
                        continuation: "stops" (< 3 mm), "runs_on" (continues without a junction),
                        "branches" (its continuation contains a junction with >= 3 arms of >= 3 mm),
                        "midline" (reaches the mid-sagittal zone, toward the sagittal sinus / falx, from an
                        end that is not itself there). Two "stops" is the cSS pattern; a vein runs on.
Tracking stays in the pial band (cortex + subarachnoid CSF): the target is veins running ALONG the surface.

A cSS focus touching an unrelated vein inherits that vein's network: read these with the images."""
import numpy as np
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra, connected_components
from skimage.morphology import skeletonize
from css_common import radius_profile

CONE_DEG = 60.0          # forward cone at each end (veins curve; cSS has no forward continuation at all)
TANGENT_MM = 3.0         # end tangent from the point this far back along the candidate
TERMINAL_MM = 3.0        # continuation below this = the candidate ends there
ARM_MIN_MM = 3.0         # a junction counts when >= 3 arms of at least this length meet
GAP_COST = 1.5           # bridging a 1-voxel gap costs 1.5 x its length

# neighbour offsets up to 2 voxels apart (2 = bridges a 1-voxel gap); forward half only
_OFF = np.array([(a, b, c) for a in range(-2, 3) for b in range(-2, 3) for c in range(-2, 3)
                 if (a, b, c) > (0, 0, 0) and max(abs(a), abs(b), abs(c)) <= 2])


def _graph(idx, shape, vox):
    """sparse graph over voxel coordinates idx (n x 3): 26-neighbours, plus gap-bridging edges"""
    lin = np.ravel_multi_index(idx.T, shape); order = np.argsort(lin); sl = lin[order]
    rows, cols, w = [], [], []
    for o in _OFF:
        q = idx + o
        ok = np.all((q >= 0) & (q < np.array(shape)), axis=1)
        ql = np.full(len(idx), -1, np.int64); ql[ok] = np.ravel_multi_index(q[ok].T, shape)
        pos = np.clip(np.searchsorted(sl, ql), 0, len(sl) - 1)
        hit = ok & (sl[pos] == ql)
        d = float(np.linalg.norm(o * vox)) * (GAP_COST if np.abs(o).max() > 1 else 1.0)
        rows.append(np.nonzero(hit)[0]); cols.append(order[pos[hit]]); w.append(np.full(hit.sum(), d))
    r, c, ww = np.concatenate(rows), np.concatenate(cols), np.concatenate(w)
    n = len(idx)
    return coo_matrix((np.r_[ww, ww], (np.r_[r, c], np.r_[c, r])), shape=(n, n)).tocsr()


def _skeleton_stats(net, vox):
    """skeleton length (mm) and junctions with >= 3 arms of >= ARM_MIN_MM"""
    if net.sum() < 3:
        return 0.0, 0
    sk = skeletonize(net) > 0
    step = float(np.mean(vox))
    length = float(sk.sum()) * step
    k = np.ones((3, 3, 3), int); k[1, 1, 1] = 0
    nb = ndi.convolve(sk.astype(int), k, mode="constant")
    junc = sk & (nb >= 3)
    if not junc.any():
        return length, 0
    jl, nj = ndi.label(junc, structure=np.ones((3, 3, 3), bool))
    arms, na = ndi.label(sk & ~ndi.binary_dilation(junc, structure=np.ones((3, 3, 3), bool)),
                         structure=np.ones((3, 3, 3), bool))
    arm_len = ndi.sum(np.ones_like(arms), arms, index=np.arange(1, na + 1)) * step if na else np.array([])
    long_arm = np.r_[False, arm_len >= ARM_MIN_MM]
    n_valid = 0
    for j, sl in enumerate(ndi.find_objects(jl), 1):
        if sl is None: continue
        big = tuple(slice(max(s.start - 2, 0), s.stop + 2) for s in sl)
        ring = ndi.binary_dilation(jl[big] == j, structure=np.ones((3, 3, 3), bool), iterations=2)
        touching = np.unique(arms[big][ring]); touching = touching[touching > 0]
        n_valid += int(long_arm[touching].sum() >= 3)
    return length, n_valid


def _end_type(cont_mm, reached, end_idx, shape, vox, midline):
    """stops / runs_on / branches / midline for one end (see module docstring)"""
    if cont_mm < TERMINAL_MM or not len(reached):
        return "stops"
    if midline is not None and not midline[tuple(end_idx)] and midline[tuple(reached.T)].any():
        return "midline"
    m = np.zeros(shape, bool); m[tuple(reached.T)] = True; m[tuple(end_idx)] = True
    return "branches" if _skeleton_stats(m, vox)[1] > 0 else "runs_on"


def local_tracking(cm, ves, vox, reach_mm, band=None, dark=None, midline=None):
    """cm: candidate mask, ves: vessel mask (same box, candidate +- reach_mm), vox: voxel size (mm),
    band: where tracking may go (detect_css.py passes the pial band: cortex + CSF, not deep white matter -
    continuation INTO the WM is vessel_wm_mm's job, and WM noise otherwise fakes continuation),
    dark: dark voxels (z < -2); continuation may also run through them, so that a vein draining into a
    dural sinus continues into it (a sinus is too wide to be a Hessian 'tube' at sigma <= 1.6 mm),
    midline: mid-sagittal zone (same box), for end_type midline"""
    path_ok = ves | dark if dark is not None else ves.copy()
    if band is not None:
        ves = ves & band; path_ok &= band
    vox = np.asarray(vox, float)
    nan = dict(cand_geodesic_mm=np.nan, local_cont_max_mm=np.nan, local_cont_min_mm=np.nan,
               local_cont_total_mm=np.nan, terminal_ends=np.nan, local_net_mm=np.nan,
               local_branch_points=np.nan, local_branch_per10mm=np.nan, end_type_max="", end_type_min="",
               net_radius_med_mm=np.nan, net_radius_p95_mm=np.nan)
    ci = np.argwhere(cm)
    if len(ci) < 3:
        return nan
    # ---- candidate geodesic diameter and its two ends (largest connected piece)
    Gc = _graph(ci, cm.shape, vox)
    d0 = dijkstra(Gc, indices=0)
    comp = np.isfinite(d0)
    if comp.sum() < 0.5 * len(ci):                      # start from a node of the biggest piece
        _, lab = connected_components(Gc)
        start = int(np.argmax(lab == np.bincount(lab).argmax())); d0 = dijkstra(Gc, indices=start)
    a = int(np.nanargmax(np.where(np.isfinite(d0), d0, -1)))
    da, pred = dijkstra(Gc, indices=a, return_predecessors=True)
    b = int(np.nanargmax(np.where(np.isfinite(da), da, -1)))
    L = float(da[b]); out = dict(nan, cand_geodesic_mm=round(L, 1))
    path = [b]
    while path[-1] != a and pred[path[-1]] >= 0:
        path.append(int(pred[path[-1]]))
    path = np.array(path[::-1])                         # a ... b along the candidate
    dpath = da[path]
    # ---- forward continuation from each end through vessel voxels only (candidate voxels excluded)
    nodes_v = np.argwhere(path_ok & ~cm)
    cont, etype = [], []
    m_ = min(TANGENT_MM, L / 2)
    back_a = path[min(np.searchsorted(dpath, m_), len(path) - 1)]
    back_b = path[max(np.searchsorted(dpath, L - m_, side="right") - 1, 0)]
    for end, back in ((a, back_a), (b, back_b)):
        pe = ci[end] * vox; t = pe - ci[back] * vox
        if np.linalg.norm(t) < 1e-6 or not len(nodes_v):
            cont.append(0.0); etype.append("stops"); continue
        t /= np.linalg.norm(t)
        rel = nodes_v * vox - pe; dist = np.linalg.norm(rel, axis=1)
        fwd = (dist <= reach_mm) & ((rel @ t >= np.cos(np.radians(CONE_DEG)) * dist) | (dist <= 2.0 * vox.max()))
        if not fwd.any():
            cont.append(0.0); etype.append("stops"); continue
        idx = np.vstack([ci[end][None], nodes_v[fwd]])
        dd = dijkstra(_graph(idx, cm.shape, vox), indices=0, limit=reach_mm)
        c_ = float(np.max(dd[np.isfinite(dd)])); cont.append(c_)
        etype.append(_end_type(c_, idx[1:][np.isfinite(dd[1:])], ci[end], cm.shape, vox, midline))
    order = np.argsort(cont)[::-1]
    out.update(end_type_max=etype[order[0]], end_type_min=etype[order[1]])
    out.update(local_cont_max_mm=round(max(cont), 1), local_cont_min_mm=round(min(cont), 1),
               local_cont_total_mm=round(sum(cont), 1),
               terminal_ends=int(sum(c < TERMINAL_MM for c in cont)))
    # ---- vessel network connected to the candidate within reach: skeleton length and junctions
    allv = np.argwhere(ves | cm)
    if len(allv) >= 3:
        G = _graph(allv, cm.shape, vox)
        src = np.nonzero(cm[tuple(allv.T)])[0]
        dn = dijkstra(G, indices=src, limit=reach_mm, min_only=True)
        net = np.zeros(cm.shape, bool); net[tuple(allv[np.isfinite(dn)].T)] = True
        length, nb = _skeleton_stats(net | cm, vox)
        out.update(local_net_mm=round(length, 1), local_branch_points=int(nb),
                   local_branch_per10mm=round(10.0 * nb / max(length, 1.0), 2))
        rn = radius_profile(path_ok | ves, vox, sample=net & ~cm)
        if len(rn):
            out.update(net_radius_med_mm=round(float(np.median(rn)), 2),
                       net_radius_p95_mm=round(float(np.percentile(rn, 95)), 2))
    return out
```

## scripts/detect_css.py

```python
#!/usr/bin/env python3
"""cSS candidate detector v4 (v4.17) - definition-driven.
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
  stays on the surface (v4.11)    vessel_wm_mm (tube continues   ~0             medullary / transcortical
                                  into the WM), vessel_ext_mm                   vein: several mm
  runs along the cortex (v4.11)   axis_normal (long axis . normal) low          high (penetrating vessel,
                                                                                vein down a sulcus)
  remote from ICH                 ich_dist_mm, near_ich_suggest  far            <=5 mm
  supratentorial, not basal       infra_frac, basal_frac,        ~0             cerebellum/tentorium,
                                  rel_height                                    basal cisterns, skull base
  supratentorial                  infratentorial                 0              1 (classical SS)
  chronic, not acute cSAH         flair_csf_z, flair_bright      FLAIR not bright (hint only:
                                                                 chronic cSS can show mild FLAIR signal)
  mimic: cortical vein thrombosis flair_ctx_z, flair_ctx_bright  FLAIR-bright cortex nearby -> check
  not a microbleed (AJNR 2016)    parenchyma_frac, extent_mm,    surface/CSF     <=10 mm and >=half
                                  cmb_like                                       in parenchyma
                                  (only if work/ID_flair_swispace.nii.gz exists)

Candidate generation (unchanged idea from v3, tightened):
  - search zone = pial band: cortex + subarachnoid CSF, -3.5..+4 mm from the pial surface, minus
    white matter / deep nuclei (v3 dilated cortex into WM, which is darker than GM on SWI)
  - darkness z vs a LOCAL cortex reference that now EXCLUDES locally dark voxels (2-pass), so
    extensive cSS no longer darkens its own reference
  - 2-D Sato ridge + hysteresis, per hemisphere; >=25 mm3, >=2 slices, elongation >=2.0
Labels: work/ID_seg_swispace.nii.gz (from T1 via prep_anat.sh when available, else SWI SynthSeg).
Lobar ICH mask: work/ID_ich.nii.gz if you drew one (freeview), else found automatically.
v4.14 bookkeeping (no change to which candidates are shown or excluded, except B1 below):
  review/ID_qc.json        case QC (css_qc.py): protocol, orientation, coverage, anatomy, pial-edge agreement,
                           status ok / low_confidence_review / insufficient_quality (reported only)
  work/ID_dropped.csv/.nii.gz  components removed by the generation gates (size, slices, elongation, no
                           cortex nearby) - before v4.14 they vanished without a record
  new columns (appended): candidates qc_status, rule_version, n_pieces; excluded also exclusion_confidence
  B1 fix: extent_mm was underestimated for candidates touching the first array slice/row/column
v4.17 (audit Stage 2): sulcal CSF = brain voxels outside the labelled tissue (recon-all aparc+aseg leaves sulcal CSF
  as label 0; label 24 alone missed it -> surface_contact and the FLAIR CSF reference were wrong on that route);
  3-D Hessian on the mask-filled image (zeros outside the mask made edge curvature); E8 artifact-zone rule needs
  T1 anatomy like E2/E3; a FLAIR older than the SWI is ignored; work/ID_gen.nii.gz records per voxel which
  generation stage it passed (label_missed.py uses it); review/ID_manifest.json records inputs (sha256),
  versions and environment.
v4.16: bounded local vein tracking columns (css_vessel.py; recorded only - no rule or score uses them):
  cand_geodesic_mm, local_cont_max/min/total_mm, terminal_ends, local_net_mm, local_branch_points/_per10mm
"""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, json, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
from skimage.filters import sato
from skimage.measure import regionprops
from skimage.morphology import skeletonize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, load_lut, fresh, sha256, update_manifest, environment_record, skeleton_length, radius_profile
import css_qc, css_vessel

subj   = sys.argv[1]
z_thr  = float(sys.argv[2]) if len(sys.argv) > 2 else -2.5
r_pct  = float(sys.argv[3]) if len(sys.argv) > 3 else 85
RIM_MM = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
RANK   = os.environ.get("CSS_RANK", "v4")
RULESET = "v4.18"                        # written to every row as rule_version
Z_LOW, R_PCT_LOW = -2.0, 80           # hysteresis (weak) thresholds
GROW_Z = -1.5                          # 'dark' map used for region growing at scoring
MIN_MM3, MIN_SLICES, MIN_ELONG = 25.0, 2, 2.0   # v4: 2.0 (tram-track = two parallel plates, ~2.4;
                                                # round microbleeds stay ~1-1.5)
PIAL_IN_MM, PIAL_OUT_MM = 3.5, 4.0     # pial band: cortex thickness + blooming / subarachnoid CSF
# v4.8 false-positive reduction (reader feedback on P006: duplicates, specks, obvious veins, skull-base
# artifact, faint normal cortex). v4.9.2: microbleed-like is a flag only; skull-base / basal rules need T1
# anatomy and a scan that includes the skull base. Excluded candidates are NOT deleted: they are listed with the reason
# in review/ID_excluded.csv (+ .nii.gz). CSS_KEEP_ALL=1 switches all exclusions off.
# v4.18: 2-D ridge (Sato) scales in mm, not pixels. [1, 2, 3] px were 0.47-1.4 mm on P010 (0.47 mm) but
# 0.86-2.6 mm on P011 (0.86 mm); these values reproduce the validated P006 behaviour (0.573 mm -> 1.05/2.09/3.14 px)
SATO_SIGMAS_MM = (0.6, 1.2, 1.8)
MERGE_GAP_MM = 3.0      # fragments of one dark line closer than this are merged into one candidate
MERGE_MAX_MM = 40.0     # ...but only while the merged object stays focus-sized: on P006 unlimited merging
                        # chained veins along a whole hemisphere into 10-17 cm objects (v4.9.1)
MIN_EXTENT_MM = 6.0     # specks: smaller than this overall are not curvilinear cSS
MIN_EDGE_SD = 0.75      # "well-defined": candidate must be darker than its 1-2 mm surroundings by this
KEEP_ALL = os.environ.get("CSS_KEEP_ALL", "") == "1"
HESS_SIGMAS_MM = (0.8, 1.6)            # 3-D shape analysis scales
# v4.11 vessel exclusion: a dark line that continues as a TUBE beyond the cortical band (into the white
# matter = medullary / transcortical vein) and runs perpendicular to the cortex is a vessel, not cSS
VES_WM_MM = 15.0        # how deep under the cortical band the vessel map reaches
MIDLINE_END_MM = 4.0    # v4.18 end_type "midline": continuation reaches within this of the mid-sagittal plane
VES_REACH_MM = 20.0     # vessel continuation is measured up to this distance (bounded: vein_tree saturated)
OFF_CTX_FRAC, OFF_PIAL_MM = 0.30, 0.30   # v4.12 "off the cortex" exclusion (T1 anatomy only)
VES_WM_Z = -3.0         # WM vessel voxels: this many SD darker than the white matter
VES_MIN_MM = 3.0        # vessel pieces shorter than this are noise
VES_CONE_DEG = 30.0     # WM continuation counts only inside this cone, straight on from the candidate
VES_EXCL_WM_MM = 4.0    # exclusion: tube continues >= this far into the WM ...
VES_EXCL_AXIS = 0.5     # ... AND the candidate runs perpendicular to the cortex (|long axis . normal|)
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
BASAL = [16, 28, 60, 17, 18, 53, 54]   # brainstem, ventral DC, hippocampus, amygdala: basal cisterns
INFRA_MM, BASAL_MM = 5.0, 5.0          # "near" = within this distance of those structures
CSF = [24]

base = _base()
swi = nib.load(f"{base}/data/{subj}_swi.nii")
I   = swi.get_fdata().astype(np.float32)
seg = nib.load(f"{base}/work/{subj}_seg_swispace.nii.gz").get_fdata().astype(np.int32)
vox = tuple(float(v) for v in swi.header.get_zooms()[:3])
# v4.9.2 anatomy confidence: millimetre position rules (skull base, basal cisterns) are only trusted
# with T1-derived labels; SynthSeg on the SWI is too coarse (P010: real-looking lesions excluded)
_t1p = f"{base}/work/{subj}_t1seg_swispace.nii.gz"
# v4.18: the slice-wise ridge filter, slice counts, merge gap and "long structure" assume array axis 2 is the
# superior-inferior / slice axis. Stop instead of running in the wrong plane (CSS_ALLOW_ORIENT=1 overrides).
_si = int(np.argmax(np.abs(swi.affine[2, :3])))
if (_si != 2 or vox[2] < max(vox[0], vox[1]) - 1e-3) and os.environ.get("CSS_ALLOW_ORIENT") != "1":
    os.makedirs(f"{base}/review", exist_ok=True)
    _why = (f"array axis 2 is not the slice / S-I axis (orientation {''.join(nib.aff2axcodes(swi.affine))}, "
            f"voxel {[round(v, 3) for v in vox]})")
    json.dump(dict(subject=subj, status="insufficient_quality", orientation_ok=False, si_axis=_si,
                   orientation="".join(nib.aff2axcodes(swi.affine)), voxel_mm=[round(v, 3) for v in vox],
                   reasons=[_why], stopped=True), open(f"{base}/review/{subj}_qc.json", "w"), indent=2)
    print(f"{subj}: STOPPED - {_why}.\n  Reslice the SWI so the slices are axial (e.g. mri_convert --out_orientation LAS"
          f" after checking the acquisition plane), or set CSS_ALLOW_ORIENT=1 to run anyway (geometry unreliable)."
          f"\n  Earlier outputs in review/ for {subj}, if any, are from a previous run.")
    sys.exit(2)
T1_ANAT = os.path.exists(_t1p) and os.path.getmtime(_t1p) >= os.path.getmtime(f"{base}/data/{subj}_swi.nii")
MERGE_GAP_MM = min(max(MERGE_GAP_MM, 2.0 * vox[2]), 4.0)   # thick slices: inter-slice gaps must not split a line
s2d = np.ones((3, 3, 1), bool); s3 = np.ones((3, 3, 3), bool)
mm2vox = lambda mm: tuple(mm / v for v in vox)

brain  = ndi.binary_fill_holes(I > 0)
cortex = (seg >= 1000) | (seg == 3) | (seg == 42)
it = max(1, int(round(2.0 / vox[0])))
edge_band = brain & ~ndi.binary_erosion(brain, structure=s2d, iterations=it)
# ---- 0. signed distance to the pial surface (mm; <0 inside tissue, >0 in CSF) + normals
tissue = ndi.binary_fill_holes(brain & (seg > 0) & ~np.isin(seg, CSF))    # fills ventricles
# v4.17: CSF = brain voxels that are not labelled tissue. SynthSeg labels extracerebral/sulcal CSF 24, but
# recon-all aparc+aseg leaves sulcal CSF unlabelled (0) and uses 24 for a small basal region only
csf_space = brain & ~tissue
surface = (~brain) | csf_space
near_surface = ndi.binary_dilation(surface, structure=s2d, iterations=it)
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
# v4.19.2 regional + contralateral reference (recorded only). LEMON P201 (healthy): 5 of 6 candidates were
# naturally dark cortex (motor strip, cuneus/precuneus: iron/myelin), patch-shaped like cSS. Physiological
# darkness covers much of its anatomical region and is roughly symmetric; cSS is focal and usually one-sided.
_ci = np.nonzero(cortex & brain); _cl = seg[_ci]; _cv = Ic[_ci].astype(np.float32); del _ci
_o = np.argsort(_cl, kind="stable"); _cl, _cv = _cl[_o], _cv[_o]
_u, _st = np.unique(_cl, return_index=True); _en = np.r_[_st[1:], len(_cl)]
REG = {}
for _c, _a, _b in zip(_u, _st, _en):
    _v = np.sort(_cv[_a:_b])
    if len(_v) >= 200:
        _m = float(np.median(_v)); REG[int(_c)] = (_v, _m, float(1.4826 * np.median(np.abs(_v - _m)) + 1e-6))
del _cl, _cv, _o
def mirror_code(c):
    if 1000 <= c < 2000: return c + 1000
    if 2000 <= c < 3000: return c - 1000
    return {3: 42, 42: 3}.get(c, -1)
def region_ref(code, vals_c, own_vals):
    """vals_c: candidate intensities; own_vals: the candidate's voxels that are in its own label's cortex"""
    nan = dict(region_z=np.nan, region_dark_frac=np.nan, mirror_region_z=np.nan, cand_mirror_z=np.nan)
    if code not in REG or not len(vals_c): return nan
    v, m, sd_ = REG[code]; cm_ = float(np.median(vals_c))
    out = dict(nan, region_z=round((cm_ - m) / sd_, 2))
    # share of the region's cortex (candidate itself left out) at least half as dark as the candidate, measured
    # from the whole-brain cortex median: high = the region is broadly dark (physiological pattern)
    thr = mu + 0.5 * (cm_ - mu)
    n_all = int(np.searchsorted(v, thr, side="right")) - int((own_vals <= thr).sum())
    n_tot = len(v) - len(own_vals)
    if cm_ < mu and n_tot > 0: out["region_dark_frac"] = round(max(n_all, 0) / n_tot, 3)
    mc = mirror_code(code)
    if mc in REG:
        vm, mm_, sdm = REG[mc]
        out["mirror_region_z"] = round((m - mm_) / np.sqrt((sd_ ** 2 + sdm ** 2) / 2), 2)
        out["cand_mirror_z"] = round((cm_ - mm_) / sdm, 2)
    return out

# ---- 2. search zone: pial band minus skull-strip rim, white matter and deep nuclei
dist_mm = ndi.distance_transform_edt(brain, sampling=vox)
rim = brain & (dist_mm < RIM_MM)
If = Ic.copy(); If[~brain] = mu
shell = brain & (pial_sd >= -PIAL_IN_MM) & (pial_sd <= PIAL_OUT_MM) & ~np.isin(seg, WM + DEEP)
shell |= edge_band & ~np.isin(seg, WM + DEEP)
zone = shell & ~rim

# ---- 3. curvilinear dark structures, slice by slice (candidate generation)
R = np.zeros_like(I)
_inpl = float(np.sqrt(vox[0] * vox[1]))
SATO_PX = [max(0.5, m / _inpl) for m in SATO_SIGMAS_MM]
if os.environ.get("CSS_SATO_LEGACY") == "1":     # v4.17 behaviour, for before/after comparisons on real data
    SATO_PX = [1.0, 2.0, 3.0]
for k in range(I.shape[2]):
    if I[:, :, k].max() > 0:
        R[:, :, k] = sato(If[:, :, k], sigmas=SATO_PX, black_ridges=True)
r_hi = np.percentile(R[shell], r_pct); r_lo = np.percentile(R[shell], R_PCT_LOW)
strong = zone & (z < z_thr) & (R > r_hi)
weak   = zone & (z < Z_LOW) & (R > r_lo)          # strong is a subset of weak
# v4.17 generation-stage record per voxel (bits): 1 search zone, 2 skull-strip rim, 4 WM / deep-nucleus label,
# 8 z < -2 (weak dark), 16 z < z_thr (strong dark), 32 ridge > weak pct, 64 ridge > strong pct,
# 128 in a component kept by hysteresis -> work/ID_gen.nii.gz, read by label_missed.py
GEN = (zone * 1 + rim * 2 + np.isin(seg, WM + DEEP) * 4 + (brain & (z < Z_LOW)) * 8 + (brain & (z < z_thr)) * 16
       + (brain & (R > r_lo)) * 32 + (brain & (R > r_hi)) * 64).astype(np.uint8)

# ---- 4. hemisphere map; hysteresis labelling separately in each hemisphere
Lm = np.isin(seg, LH) | ((seg >= 1000) & (seg < 2000))
Rm = np.isin(seg, RH) | (seg >= 2000)
hemiL = ndi.distance_transform_edt(~Lm, sampling=vox) <= ndi.distance_transform_edt(~Rm, sampling=vox)
dmid = np.where(hemiL, ndi.distance_transform_edt(hemiL, sampling=vox),
                ndi.distance_transform_edt(~hemiL, sampling=vox))
lab = np.zeros(I.shape, np.int32); n_raw = 0; nxt = 0
n_pieces = None; n_before_merge = None     # set by the merge step
for side in (hemiL, ~hemiL):
    wl, nw = ndi.label(weak & side, structure=s3)
    n_raw += nw
    keep = np.unique(wl[strong & side]); keep = keep[keep > 0]
    remap = np.zeros(nw + 1, np.int32)
    remap[keep] = np.arange(1, len(keep) + 1) + nxt
    nz = wl > 0
    lab[nz] = remap[wl[nz]]
    nxt += len(keep)
GEN[lab > 0] |= 128
nib.save(nib.Nifti1Image(GEN, swi.affine), f"{base}/work/{subj}_gen.nii.gz"); del GEN
del R, wl, weak, strong

# ---- 4b. merge fragments of the same structure: a line broken by small gaps (noise, partial volume,
# slice-to-slice wobble) was reported as 2-3 separate candidates; also joins the two banks of a
# tram-track. Done per hemisphere so left and right never merge.
if MERGE_GAP_MM > 0 and lab.max() > 0:
    near_lab = ndi.distance_transform_edt(lab == 0, sampling=vox) <= MERGE_GAP_MM / 2
    merged = np.zeros_like(lab); nxt = 0
    for side in (hemiL, ~hemiL):
        gl, ng = ndi.label(near_lab & side, structure=s3)
        sel = (lab > 0) & side
        merged[sel] = gl[sel] + nxt
        nxt += ng
    # undo merges that would chain into an object longer than MERGE_MAX_MM: keep the pieces separate
    gmax = int(merged.max()); big = np.zeros(gmax + 1, bool)
    for g, sl in enumerate(ndi.find_objects(merged), 1):
        if sl is not None and np.sqrt(sum(((s_.stop - s_.start) * v) ** 2 for s_, v in zip(sl, vox))) > MERGE_MAX_MM:
            big[g] = True
    n_chains = int(big.sum())
    merged = np.where(big[merged], lab + gmax, merged)
    u = np.unique(merged[merged > 0]); remap = np.zeros(int(merged.max()) + 1, np.int32)
    remap[u] = np.arange(1, len(u) + 1)
    n_before_merge = int(len(np.unique(lab[lab > 0])))
    _orig = lab
    lab = remap[merged]
    _sel = _orig > 0
    _pairs = np.unique(np.stack([lab[_sel], _orig[_sel]]), axis=1)
    n_pieces = np.bincount(_pairs[0], minlength=int(lab.max()) + 1)    # original pieces per candidate
    del _orig, _sel, _pairs
    del near_lab, merged

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
            h = ndi.gaussian_filter(If, sg, order=order)[idx] / (vox[a] * vox[b]) * s_mm ** 2   # v4.17: If
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
# v4.11 vessel map that ALSO covers the white matter under the cortex: veins do not stop at the cortex -
# medullary / transcortical veins run radially into the WM, cortical veins continue along the surface.
# cSS is confined to the surface. WM is darker than cortex on SWI, so WM voxels are judged against WM.
_wmm = np.isin(seg, [2, 41])
wm_mu = float(np.median(Ic[_wmm])) if _wmm.any() else mu
wneed = _wmm & ((Ic - wm_mu) / sd < VES_WM_Z) & (pial_sd > -VES_WM_MM - PIAL_IN_MM) & ~need
if wneed.any():
    widx = np.nonzero(wneed); wbest = np.full(len(widx[0]), -np.inf, np.float32)
    weig = np.zeros((len(widx[0]), 3), np.float32)
    for s_mm in HESS_SIGMAS_MM:
        sg = mm2vox(s_mm); H = np.zeros((len(widx[0]), 3, 3), np.float32)
        for a in range(3):
            for b in range(a, 3):
                order = [0, 0, 0]; order[a] += 1; order[b] += 1
                h = ndi.gaussian_filter(If, sg, order=order)[widx] / (vox[a] * vox[b]) * s_mm ** 2
                H[:, a, b] = h; H[:, b, a] = h
        w = np.linalg.eigvalsh(H); w = np.take_along_axis(w, np.argsort(np.abs(w), axis=1), 1)
        upd = w[:, 2] > wbest; wbest[upd] = w[upd, 2]; weig[upd] = w[upd]
        del H, w
    # tube = two large curvatures and one small one; a round blob (haematoma, microbleed) is not a vessel
    _l3 = np.maximum(weig[:, 2], 1e-6)
    TUBE[widx] = np.where((weig[:, 2] > 0) & (np.abs(weig[:, 0]) / _l3 < 0.5), np.abs(weig[:, 1]) / _l3, np.nan)
    del weig, wbest
VES = ((need & (z < Z_LOW)) | wneed) & (TUBE > 0.5)        # dark tubes: band + WM
del wneed
# noise specks are not vessels: keep tube pieces at least VES_MIN_MM long
_vl, _ = ndi.label(VES, structure=s3)
_ok = np.zeros(int(_vl.max()) + 1, bool)
for n, sl_ in enumerate(ndi.find_objects(_vl), 1):
    if sl_ is not None:
        _ok[n] = np.sqrt(sum(((s_.stop - s_.start) * v) ** 2 for s_, v in zip(sl_, vox))) >= VES_MIN_MM
VES = _ok[_vl]; del _vl, _ok
deep_wm = _wmm & (ndi.distance_transform_edt(_wmm, sampling=vox) >= 1.5)   # WM proper, not the cortex border
# v4.16 dark structures for local vein tracking: dark voxels in pieces >= VES_MIN_MM (isolated noise voxels
# below z -2 would otherwise chain up through the tracker's gap bridging); includes wide dark structures
# such as a dural sinus, which are not Hessian tubes at sigma <= 1.6 mm
_dl, _ = ndi.label(brain & (z < Z_LOW), structure=s3)
_ok = np.zeros(int(_dl.max()) + 1, bool)
for n, sl_ in enumerate(ndi.find_objects(_dl), 1):
    if sl_ is not None:
        _ok[n] = np.sqrt(sum(((s_.stop - s_.start) * v) ** 2 for s_, v in zip(sl_, vox))) >= VES_MIN_MM
DARKS = _ok[_dl]; del _dl, _ok
nib.save(nib.Nifti1Image(VES.astype(np.uint8), swi.affine), f"{base}/work/{subj}_vessels.nii.gz")   # for freeview
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
brain_in = ndi.binary_erosion(brain, structure=s3)      # one voxel inside the mask edge
def mirror_dark(pts):
    p = pts * vox; q = p - 2 * ((p - pc) @ mvec)[:, None] * mvec[None, :]
    qi = np.round(q / vox).astype(int)
    ok = np.all((qi >= 0) & (qi < np.array(I.shape)), axis=1)
    qi = qi[ok]; ok[ok] = brain_in[qi[:, 0], qi[:, 1], qi[:, 2]]
    # v4.18: mirror points outside the brain read as dark (stripped background) -> ignored; NaN if most are outside
    if ok.mean() < 0.5:
        return float("nan")
    qi = np.round(q[ok] / vox).astype(int)
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
FL = None; fl_note = ""
if os.path.exists(fl_p) and not fresh(fl_p, f"{base}/data/{subj}_swi.nii"):
    fl_note = "FLAIR IGNORED (older than the SWI - re-run prep_anat.sh)"
elif os.path.exists(fl_p):
    FL = nib.load(fl_p).get_fdata().astype(np.float32)
    csf_ref = FL[csf_space & (FL > 0)]                      # v4.17 (was label 24 only)
    tis_ref = FL[tissue & brain & (FL > 0)]
    ctx_ref = FL[cortex & brain & (FL > 0)]
    if csf_ref.size and tis_ref.size and ctx_ref.size:
        f_mu = float(np.median(csf_ref)); f_sd = float(1.4826 * np.median(np.abs(tis_ref - np.median(tis_ref))) + 1e-6)
        fc_mu = float(np.median(ctx_ref))
    else:
        FL = None

# ---- 10. features per candidate
cortex_dist = ndi.distance_transform_edt(~cortex, sampling=vox)      # 0 inside cortex
# v4.9 normal anatomy: cerebellum/brainstem (tentorial interface) and basal cisterns
d_infra = ndi.distance_transform_edt(~np.isin(seg, INFRA), sampling=vox) if np.isin(seg, INFRA).any() \
    else np.full(I.shape, np.inf, np.float32)
d_basal = ndi.distance_transform_edt(~np.isin(seg, BASAL), sampling=vox) if np.isin(seg, BASAL).any() \
    else np.full(I.shape, np.inf, np.float32)
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
        bp += int((sk & (nb >= 3)).sum()); length += skeleton_length(sk, vox)
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
# v4.14: components removed by the generation gates are logged (work/ID_dropped.csv/.nii.gz), not shown
dropped = []
def drop(p, why, vol, n_sl=None, elong=None):
    c_ = np.asarray(p.centroid) / np.array(vox)              # regionprops centroid is in mm (spacing)
    dropped.append(dict(drop_id=len(dropped) + 1, reason=why, volume_mm3=round(vol, 1),
                        n_slices=n_sl, elongation=None if elong is None else round(elong, 1),
                        vox_i=int(round(c_[0])), vox_j=int(round(c_[1])), vox_k=int(round(c_[2])),
                        hemi="L" if hemiL[tuple(np.clip(np.round(c_).astype(int), 0, np.array(I.shape) - 1))] else "R",
                        _lab=p.label))
for p in regionprops(lab, spacing=vox):
    vol = float(p.area)                         # already mm3 because spacing is given
    if vol < MIN_MM3:
        drop(p, f"volume <{MIN_MM3:g} mm3", vol); continue
    sl = objs[p.label - 1]
    pad = tuple(slice(max(s.start - 3, 0), s.stop + 3) for s in sl)
    mm = lab[pad] == p.label
    n_sl = int(mm.any(axis=(0, 1)).sum())
    if n_sl < MIN_SLICES:
        drop(p, f"<{MIN_SLICES} slices", vol, n_sl); continue
    elong = p.axis_major_length / max(p.axis_minor_length, 1e-3)
    # elongation rejects round microbleeds, which are single pieces; a merged candidate (two banks,
    # a bent line) can be compact overall and is not rejected for it (v4.9.2)
    if elong < MIN_ELONG and (n_pieces is None or n_pieces[p.label] <= 1):
        drop(p, f"elongation <{MIN_ELONG:g}", vol, n_sl, elong); continue
    seg_nb = seg[pad][ndi.binary_dilation(mm, iterations=2)]
    nb = seg_nb[seg_nb >= 1000]
    # infratentorial only when clearly cerebellar/brainstem with little cerebral cortex nearby, so
    # occipital / inferior temporal cSS next to the tentorium is not dropped from the 0-4 score
    infra = bool(np.isin(seg_nb, INFRA).mean() > 0.5 and nb.size < 0.2 * seg_nb.size)
    if nb.size == 0 and not infra:
        drop(p, "no cortex label within 2 voxels", vol, n_sl, elong); continue
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
    # long axis vs the cortical surface normal: cSS runs ALONG the cortex (~0), penetrating / bridging
    # vessels and veins running down the middle of a sulcus run ALONG the normal (~1)
    pm = pts * np.array(vox); pm = pm - pm.mean(0)
    ua = np.linalg.svd(pm[:: max(1, len(pm) // 2000)], full_matrices=False)[2][0] if len(pm) >= 3 else np.zeros(3)
    axn = float(np.median(np.abs(sum(ua[a] * normal[a][tuple(pts.T)] for a in range(3)))))
    # v4.11 vessel continuation, bounded to VES_REACH_MM: follow dark TUBES connected to the candidate
    #  vessel_ext_mm = how far the tube continues beyond the candidate (cortical vein running on)
    #  vessel_wm_mm  = how far it runs on, in line with the candidate, inside the white matter (medullary vein)
    pv = tuple(int(np.ceil(VES_REACH_MM / v)) for v in vox)
    vb = tuple(slice(max(s_.start - q_, 0), min(s_.stop + q_, n_)) for s_, q_, n_ in zip(sl, pv, I.shape))
    cm = lab[vb] == p.label
    dcm = ndi.distance_transform_edt(~cm, sampling=vox)
    vv_ = VES[vb] & (dcm <= VES_REACH_MM)
    # 1-voxel gaps are bridged (thin veins break up at the cortex/WM contrast step and in noise);
    # distances are then measured on real vessel voxels only
    gl_, _ = ndi.label(ndi.binary_dilation(vv_, structure=s3) | cm, structure=s3)
    gid = np.unique(gl_[cm]); gid = gid[gid > 0]
    reach = np.isin(gl_, gid) & vv_ & ~cm
    ves_ext = float(dcm[reach].max()) if reach.any() else 0.0
    rwm = reach & deep_wm[vb] & (ich_dist[vb] > 3.0)          # a haematoma's rim is not a vessel
    # only WM vessel that carries STRAIGHT ON beyond the candidate's deep end:
    # a transcortical vein continues in line; cSS that merely touches an unrelated vessel does not
    # cone (half-angle VES_CONE_DEG) from the candidate's deepest point, pointing inward: the mean
    # inward surface normal + the centre->deepest-point direction (where a transcortical vein exits)
    if rwm.any():
        vx_ = np.array(vox); deep = pts[np.argmin(psd)] * vx_
        nin = -np.array([normal[a][tuple(pts.T)].mean() for a in range(3)])
        dvec = deep - pts.mean(0) * vx_
        ud = nin / (np.linalg.norm(nin) + 1e-6) + dvec / (np.linalg.norm(dvec) + 1e-6)
        ud = ud / np.linalg.norm(ud) if np.linalg.norm(ud) > 0.3 else ua
        fp = (np.argwhere(rwm) + [s_.start for s_ in vb]) * vx_ - deep
        inline = (fp @ ud) >= np.cos(np.radians(VES_CONE_DEG)) * np.linalg.norm(fp, axis=1)
        ves_wm = float(dcm[rwm][inline].max()) if inline.any() else 0.0
    else:
        ves_wm = 0.0
    # vessel_run_mm (recorded only): tube continuing straight on beyond EITHER end of the candidate's long
    # axis (a cortical vein runs on towards a sinus; cSS stops). Bounded by VES_REACH_MM.
    ves_run = 0.0
    if reach.any() and len(pts) >= 3:
        vx_ = np.array(vox); pm_ = pts * vx_; proj = (pm_ - pm_.mean(0)) @ ua
        fq = (np.argwhere(reach) + [s_.start for s_ in vb]) * vx_
        for end, sgn in ((pm_[np.argmax(proj)], 1.0), (pm_[np.argmin(proj)], -1.0)):
            fe = fq - end; ok_ = (fe @ (sgn * ua)) >= np.cos(np.radians(VES_CONE_DEG)) * np.linalg.norm(fe, axis=1)
            if ok_.any(): ves_run = max(ves_run, float(np.linalg.norm(fe[ok_], axis=1).max()))
    del dcm, gl_, reach, rwm
    # v4.16 bounded local vein tracking along the vessel map (recorded only, css_vessel.py)
    vt = css_vessel.local_tracking(cm, VES[vb], vox, VES_REACH_MM, band=(pial_sd[vb] > -PIAL_IN_MM) & ~np.isin(seg[vb], WM + DEEP),
                                   dark=DARKS[vb], midline=dmid[vb] < MIDLINE_END_MM)
    mdark = mirror_dark(pts)
    _own = mm & (seg[pad] == code) & cortex[pad]
    rref = region_ref(code, Ic[pad][mm], Ic[pad][_own])
    # v4.18 thickness (recorded only): half-width along the candidate's centre line
    rp = radius_profile(mm, vox)
    rad = dict(cand_radius_med_mm=round(float(np.median(rp)), 2) if len(rp) else np.nan,
               cand_radius_p95_mm=round(float(np.percentile(rp, 95)), 2) if len(rp) else np.nan,
               cand_radius_cv=round(float(np.std(rp) / np.mean(rp)), 2) if len(rp) and np.mean(rp) > 0 else np.nan)
    # microbleed vs cSS (Charidimou, AJNR 2016;37:E43): microbleeds are small (generally 2-5 mm),
    # round/oval and at least half surrounded by brain parenchyma - cSS lies on the surface / in CSF
    shell1 = ndi.binary_dilation(mm, structure=s3) & ~mm
    pfrac = float(tissue[pad][shell1].mean()) if shell1.any() else np.nan
    # bounding-box diagonal of the candidate itself (v4.14 B1: was pad width - 6, too small at the volume edge)
    ext_mm = float(np.sqrt(sum(((s_.stop - s_.start) * v) ** 2 for s_, v in zip(sl, vox))))
    cmb_like = bool(np.isfinite(pfrac) and pfrac >= 0.5 and ext_mm <= 10.0)
    # "well-defined" (consensus definition): darker than the 1-2 mm surroundings by >= MIN_EDGE_SD;
    # normal dark cortex fades gradually into its neighbourhood
    ring2 = ndi.binary_dilation(mm, structure=s3, iterations=2) & ~mm & brain[pad]
    edge_c = float(np.median(z[pad][ring2]) - np.median(z[pad][mm])) if ring2.any() else np.nan
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
    # mirror_dark_frac is recorded only (venous anatomy can be asymmetric, cSS can be bilateral)
    for flag, f in ((artz, 0.5), (midz, 0.7),
                    (longs, 0.5), (infra, 0.3)):
        if flag: s4 *= f
    # two separate, rule-based evidence summaries for the reader (0-1, NOT probabilities -
    # untrained; to be replaced by a model fitted on expert labels):
    #  cSS  = on the surface, sheet-like, parallel to the cortex, tram-track, on a bank
    #  vein = tube-like, out in the sulcal CSF, not parallel to the cortex
    fin = lambda v: v if np.isfinite(v) else None
    ce = [on_surf, 1 - tub if np.isfinite(tub) else None, fin(aln), fin(tram), bankf]
    ve = [fin(tub), 1 / (1 + np.exp(-(pdist - 0.3) / 0.5)), 1 - aln if np.isfinite(aln) else None,
          axn, min(ves_wm / VES_EXCL_WM_MM, 1.0)]
    css_ev = float(np.mean([v for v in ce if v is not None]))
    vein_ev = float(np.mean([v for v in ve if v is not None]))
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
                     parenchyma_frac=round(pfrac, 2), extent_mm=round(ext_mm, 1), cmb_like=int(cmb_like),
                     near_ich_suggest=int(idist <= 5.0), infratentorial=int(infra),
                     flair_csf_z=round(fz, 2), flair_bright=int(np.isfinite(fz) and fz > 3.0),
                     flair_ctx_z=round(fcz, 2), flair_ctx_bright=int(np.isfinite(fcz) and fcz > 3.0),
                     css_evidence=round(css_ev, 2), vein_evidence=round(vein_ev, 2),
                     edge_contrast=round(edge_c, 2), _wz=float((swi.affine @ np.r_[c0, 1])[2]),
                     infra_frac=round(float((d_infra[pad][mm] <= INFRA_MM).mean()), 2),
                     basal_frac=round(float((d_basal[pad][mm] <= BASAL_MM).mean()), 2),
                     vessel_ext_mm=round(ves_ext, 1), vessel_wm_mm=round(ves_wm, 1), vessel_run_mm=round(min(ves_run, VES_REACH_MM), 1), axis_normal=round(axn, 2),
                     score_v3=round(s3_, 2), score_v4=round(s4, 2), _lab=p.label,
                     n_pieces=int(n_pieces[p.label]) if n_pieces is not None else 1, **vt, **rad, **rref))

# CSV columns up to "long_structure" are consumed by export_review.py / the workbook: keep them
cols = ["cand_id", "hemi", "region", "label", "volume_mm3", "n_slices", "elongation",
        "surface_contact", "darkness_z", "artifact_zone", "midline_zone", "score",
        "vox_i", "vox_j", "vox_k", "cortex_frac", "cortex_dist_mm", "branch_per10mm", "surface_gradient",
        "vein_like", "long_structure",
        "pial_dist_mm", "bank_frac", "tube_ratio", "surface_alignment", "tram_frac", "vein_tree_mm",
        "mirror_dark_frac", "ich_dist_mm", "near_ich_suggest", "infratentorial", "flair_csf_z",
        "parenchyma_frac", "extent_mm", "cmb_like",
        "flair_bright", "flair_ctx_z", "flair_ctx_bright", "css_evidence", "vein_evidence",
        "edge_contrast", "infra_frac", "basal_frac", "rel_height",
        "vessel_ext_mm", "vessel_wm_mm", "vessel_run_mm", "axis_normal", "score_v3", "score_v4", "accept",
        "qc_status", "rule_version", "n_pieces",                                   # v4.14, appended
        "cand_geodesic_mm", "local_cont_max_mm", "local_cont_min_mm", "local_cont_total_mm",  # v4.16,
        "terminal_ends", "local_net_mm", "local_branch_points", "local_branch_per10mm",       # recorded only
        "end_type_max", "end_type_min", "review_tier", "tier_reason",
        "cand_radius_med_mm", "cand_radius_p95_mm", "cand_radius_cv", "net_radius_med_mm", "net_radius_p95_mm",                        # v4.18, recorded only
        "region_z", "region_dark_frac", "mirror_region_z", "cand_mirror_z"]                         # v4.19.2, recorded only

# ---- 11. exclusions (v4.8) - each with an explicit, definition-based reason
# height is measured within the CEREBRUM (cortex + cerebral WM), so 0 = floor of the temporal/frontal
# lobes (skull base), 1 = vertex; the cerebellum below does not stretch the scale
_cer = (cortex | np.isin(seg, [2, 41]))[::4, ::4, ::4]
_src = _cer if _cer.any() else brain[::4, ::4, ::4]
_bz = (swi.affine @ np.c_[np.argwhere(_src) * 4, np.ones(int(_src.sum()))].T)[2]
z_lo, z_hi = float(_bz.min()), float(_bz.max())
for r in rows: r["rel_height"] = round((r["_wz"] - z_lo) / max(z_hi - z_lo, 1e-6), 2)
# coverage: a slab that is cut off at the bottom has no skull base in it - "lowest 20 %" would then
# remove ordinary cortex (P010: lateral occipital lesion on slice 6)
_zax = int(np.argmax(np.abs(swi.affine[2, :3])))
_area = brain.sum(axis=tuple(a for a in range(3) if a != _zax)).astype(float)
_nz = np.nonzero(_area)[0]
_bot = _nz[0] if swi.affine[2, _zax] > 0 else _nz[-1]
FULL_BOTTOM = bool(len(_nz) and _area[_bot] < 0.2 * _area.max())
# v4.14: each rule has an id + the version that introduced or last changed it (rule_version column) and a
# confidence for the excluded list. Confidence is a reading aid, not validated: "low" marks rules with
# little real-data support or a case whose QC status is not ok (see review/ID_qc.json).
RULE_CONF = {"E1": "high" if T1_ANAT else "medium", "E2": "medium", "E3": "medium", "E4": "high",
             "E5": "medium", "E6": "medium", "E7": "medium", "E8": "low", "E9": "medium"}
def exclusion(r):
    """-> (rule_version, reason); ("", "") = kept"""
    if KEEP_ALL: return "", ""
    if r["infratentorial"] or r["infra_frac"] >= (0.3 if T1_ANAT else 0.5):
        return "E1@v4.9.2", "infratentorial / tentorial interface (cerebellum, brainstem) - not cSS by definition"
    if T1_ANAT and r["basal_frac"] >= 0.3:
        return "E2@v4.9.2", "basal cisterns (circle of Willis / basal veins next to brainstem, mesial temporal lobe)"
    if T1_ANAT and FULL_BOTTOM and r["rel_height"] < 0.2:
        return "E3@v4.9.2", "skull base (lowest 20 % of the cerebrum)"
    if r["extent_mm"] < MIN_EXTENT_MM: return "E4@v4.8", f"speck (<{MIN_EXTENT_MM:g} mm)"
    tub_, tram_ = r["tube_ratio"], r["tram_frac"]
    if r["pial_dist_mm"] > 0.5 and np.isfinite(tub_) and tub_ > 0.45 and not (np.isfinite(tram_) and tram_ >= 0.2):
        return "E5@v4.8", "vein: tubular, in the middle of the sulcal CSF"
    if r["vessel_wm_mm"] >= VES_EXCL_WM_MM and r["axis_normal"] >= VES_EXCL_AXIS \
            and not (np.isfinite(tram_) and tram_ >= 0.4):     # clear tram-track (cSS ~0.9) overrides
        return "E6@v4.11", "vessel: continues as a tube into the white matter, runs perpendicular to the cortex"
    # v4.12 off the cortex (T1 anatomy only): cSS coats the pial surface, so it always touches the cortex.
    # A dark line lying wholly out in the sulcal CSF is a vein. Known P006 cSS: cortex_frac >= 0.55,
    # pial_dist <= -0.57; P011 false positives: median 0.16 / +1.22 mm -> thresholds with a wide margin
    if T1_ANAT and r["cortex_frac"] < OFF_CTX_FRAC and r["pial_dist_mm"] > OFF_PIAL_MM \
            and not (np.isfinite(tram_) and tram_ >= 0.2):
        return "E7@v4.12", "off the cortex: lies in the sulcal CSF without touching the cortex (vein)"
    if T1_ANAT and FULL_BOTTOM and r["artifact_zone"] and r["rel_height"] < 0.4:   # v4.17: T1 only, like E2/E3
        return "E8@v4.17", "skull-base susceptibility artifact zone"
    if np.isfinite(r["edge_contrast"]) and r["edge_contrast"] < MIN_EDGE_SD:
        return "E9@v4.8", "faint / ill-defined (normal dark cortex)"
    return "", ""
# ---- 11b. case QC (v4.14, reported only - changes no candidate)
QC = css_qc.case_qc(base, subj, swi, I, If, seg, brain, tissue, cortex, pial_sd, normal, rim, hemiL, vox,
                    T1_ANAT, FULL_BOTTOM)
for r in rows:
    rv, r["_excl"] = exclusion(r)
    r["qc_status"] = QC["status"]
    r["rule_version"] = rv or RULESET
    r["exclusion_confidence"] = ("low" if QC["status"] != "ok" else RULE_CONF[rv.split("@")[0]]) if rv else ""
    # v4.18 review tier (audit S3), RECORDED ONLY: 2 = the candidate's own flags point to a mimic. Nothing is hidden
    # and the order is unchanged (review_css.py --by-tier is opt-in). Case QC is not used: it applies to every
    # candidate of the case alike. Not validated - check the known lesions' tiers before ordering by it.
    why = [n for n, f in (("artifact zone", r["artifact_zone"]), ("vein-like", r["vein_like"]),
                          ("microbleed-like", r["cmb_like"]),
                          ("vein evidence > cSS evidence", r["vein_evidence"] > r["css_evidence"])) if f]
    r["review_tier"] = 2 if why else 1; r["tier_reason"] = "; ".join(why)
excluded = [r for r in rows if r["_excl"]]
rows = [r for r in rows if not r["_excl"]]
os.makedirs(f"{base}/review", exist_ok=True)
xl = np.zeros(I.shape, np.int16)
if excluded:
    xdf = pd.DataFrame(excluded).sort_values("score", ascending=False).reset_index(drop=True)
    xdf.insert(0, "excl_id", np.arange(1, len(xdf) + 1))
    xr = np.zeros(lab.max() + 1, np.int16); xr[xdf["_lab"].values] = xdf["excl_id"].values; xl = xr[lab]
    xdf.rename(columns={"_excl": "reason"})[["excl_id", "reason"] + [c for c in cols if c not in ("cand_id", "accept")]
                                            + ["exclusion_confidence"]] \
        .to_csv(f"{base}/review/{subj}_excluded.csv", index=False)
else:
    pd.DataFrame(columns=["excl_id", "reason"]).to_csv(f"{base}/review/{subj}_excluded.csv", index=False)
nib.save(nib.Nifti1Image(xl, swi.affine), f"{base}/review/{subj}_excluded.nii.gz")
# generation-gate drops (v4.14): not for review, but lets stress_test / check_known say why a lesion is missing
dl = np.zeros(I.shape, np.int16)
if dropped:
    ddf = pd.DataFrame(dropped)
    dr = np.zeros(lab.max() + 1, np.int16); dr[ddf["_lab"].values] = ddf["drop_id"].values; dl = dr[lab]
    ddf.drop(columns="_lab").to_csv(f"{base}/work/{subj}_dropped.csv", index=False)
else:
    pd.DataFrame(columns=["drop_id", "reason"]).to_csv(f"{base}/work/{subj}_dropped.csv", index=False)
nib.save(nib.Nifti1Image(dl, swi.affine), f"{base}/work/{subj}_dropped.nii.gz")
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
print(f"{subj}: anatomy {'T1' if T1_ANAT else 'SWI-only (position rules relaxed)'}, "
      f"coverage {'includes skull base' if FULL_BOTTOM else 'cut off at the bottom (skull-base rules off)'}, "
      f"merge gap {MERGE_GAP_MM:.1f} mm")
print(f"{subj}: {n_raw} raw components -> {len(df)} candidates "
      f"[{len(excluded)} excluded: " + (", ".join(f"{v} {k.split(' (')[0].split(':')[0]}" for k, v in
      pd.Series([r['_excl'] for r in excluded]).value_counts().items()) if excluded else "none")
      + (f"; {n_before_merge - int(lab.max())} pieces joined, {n_chains} long chains kept apart"
         if n_before_merge is not None else "") + "] "
      f"({cnt('artifact_zone')} artifact zone, {cnt('midline_zone')} midline, {cnt('vein_like')} vein-like, "
      f"{cnt('near_ich_suggest')} near ICH [{ich_src}], {cnt('infratentorial')} infratentorial)  rank={RANK}"
      + ("  FLAIR used" if FL is not None else ""))
if fl_note: print(f"{subj}: {fl_note}")
print(f"{subj}: QC {QC['status']}" + (f" ({'; '.join(QC['reasons'])})" if QC["reasons"] else "")
      + f", label shift {QC['label_shift_mm']} mm, edge bias {QC['edge_bias_mm']} mm, scatter {QC['edge_scatter_mm']} mm, {len(dropped)} components below the generation gates"
      + " [review/" + subj + "_qc.json; reported only]")
print(f"{subj}: review tier (recorded only) 1: {sum(r['review_tier'] == 1 for r in rows)}, "
      f"2: {sum(r['review_tier'] == 2 for r in rows)}\n")
if len(df):
    print(df.head(15)[["cand_id", "hemi", "region", "volume_mm3", "darkness_z", "pial_dist_mm",
                       "tube_ratio", "surface_alignment", "tram_frac", "vein_tree_mm", "score",
                       "vox_i", "vox_j", "vox_k"]].to_string(index=False))
# v4.17 provenance: inputs (sha256), parameters, versions -> review/ID_manifest.json
_inputs = {k: p_ for k, p_ in (("swi", f"{base}/data/{subj}_swi.nii"), ("labels", f"{base}/work/{subj}_seg_swispace.nii.gz"),
                               ("t1_labels", _t1p), ("destrieux", f"{base}/work/{subj}_a2009s_swispace.nii.gz"),
                               ("flair", fl_p), ("phase", f"{base}/data/{subj}_phase.nii.gz"),
                               ("registration", f"{base}/work/{subj}_swi2t1.lta"), ("ich_drawn", ich_p))
           if os.path.exists(p_)}
update_manifest(base, subj, "detect", dict(
    pipeline_version=RULESET, run_at=pd.Timestamp.now().isoformat(timespec="seconds"), args=sys.argv[1:],
    parameters=dict(z_thr=z_thr, ridge_pct=r_pct, rim_mm=RIM_MM, rank=RANK, keep_all=KEEP_ALL,
                    sato_sigmas_mm=list(SATO_SIGMAS_MM), sato_sigmas_px=[round(v, 2) for v in SATO_PX]),
    anatomy="T1" if T1_ANAT else "SWI-only", flair_used=FL is not None, qc_status=QC["status"],
    n_candidates=int(len(df)), n_excluded=len(excluded), n_dropped=len(dropped),
    inputs={k: dict(path=p_, sha256=sha256(p_), mtime=os.path.getmtime(p_)) for k, p_ in _inputs.items()},
    environment=environment_record()))
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

## scripts/expert_import.py

```python
#!/usr/bin/env python3
"""Import an expert's answers from the blinded sheet (expert_sheet.py) as REFERENCE calls.
usage: expert_import.py SUBJECT READER --css 3,17,40 [--vein 1,2] [--ich 9] [--unsure 5]
       expert_import.py SUBJECT READER --letters "1V 2VS:4 3C 4N 5U ..."   (codes from the grid)
       options: --session N (default 1; a later blinded re-read of the same sheet = 2, ...)  --score
Codes (css_common.LABELS): C, V / VS / VC / VT / VX, TV, A / AM, N / NI, MB, CA, SAH, LN, HI, IS, H, U,
optionally followed by :1-:5 confidence. Sheet numbers not listed are recorded as N (defaulted=1).
Writes
  review/ID_expert_READER.csv   candidates: cand_id, call (coarse: cSS / Vein / Normal / Artifact /
                                Near ICH / Unsure), fingerprint, then label, confidence, defaulted (v4.15)
  labels/candidate_labels.csv   long table, one row per study x reader x session x object, INCLUDING
                                excluded objects when the sheet mixed them in (expert_sheet --include-excluded)
The reader's own calls in review/ID_calls.csv are left untouched. Then scores the expert calls:
  expert_import.py ... --score     (runs mark_css.py + score_css.py with the expert's cSS / ICH calls)"""
import sys, os, re, datetime, subprocess, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, FP, match_calls, LABELS, parse_answer, cand_uid, append_rows

if len(sys.argv) < 3: sys.exit(__doc__)
subj, reader = sys.argv[1], sys.argv[2]
base = _base()
key = pd.read_csv(f"{base}/review/{subj}_expert_key.csv")
if "list" not in key.columns:                     # key written before v4.15: candidates only
    key["list"], key["obj_id"] = "candidate", key["cand_id"]
arg = lambda f: sys.argv[sys.argv.index(f) + 1] if f in sys.argv else ""
nums = lambda s: [int(x) for x in re.split(r"[,\s]+", s.strip()) if x]
session = int(arg("--session") or 1)
ans = {n: ("N", None, 1) for n in key.sheet_no}  # sheet_no -> (code, confidence, defaulted)
if arg("--letters"):
    for tok in arg("--letters").split():
        try:
            n, code, conf = parse_answer(tok)
        except ValueError as e:
            sys.exit(str(e))
        ans[n] = (code, conf, 0)
for flag, code in (("--css", "C"), ("--vein", "V"), ("--ich", "H"), ("--unsure", "U")):
    for n in nums(arg(flag)): ans[n] = (code, None, 0)
bad = [n for n in ans if n not in set(key.sheet_no)]
if bad: sys.exit(f"sheet numbers not on the sheet: {bad}")

# re-attach each sheet number to the CURRENT detector output by fingerprint (the detector may have been
# re-run since the sheet was made); the 'call' passed to match_calls is the sheet number itself
mapped, dropped = {}, 0
for lst, path, idc in (("candidate", f"{base}/review/{subj}_candidates.csv", "cand_id"),
                       ("excluded", f"{base}/review/{subj}_excluded.csv", "excl_id")):
    k_ = key[key.list == lst]
    if not len(k_): continue
    cur = pd.read_csv(path) if os.path.exists(path) else pd.DataFrame(columns=[idc] + FP)
    m, d = match_calls(cur.rename(columns={idc: "cand_id"}), k_[FP].assign(call=k_.sheet_no.values))
    dropped += d
    for oid, n in (m or {}).items():
        mapped[(lst, oid)] = (int(n), cur.set_index(idc).loc[oid])
if dropped: print(f"WARNING: {dropped} answers match no current object (detector re-run since the sheet)")

ver = "v4.13 or earlier"
cand = pd.read_csv(f"{base}/review/{subj}_candidates.csv")
if "rule_version" in cand.columns and len(cand): ver = str(cand.rule_version.iloc[0])
now = datetime.datetime.now().isoformat(timespec="seconds")
rows, long_rows = [], []
for (lst, oid), (n, r) in sorted(mapped.items(), key=lambda kv: kv[1][0]):
    code, conf, dflt = ans[n]; label, call = LABELS[code]
    if lst == "candidate":
        rows.append(dict(cand_id=oid, call=call, **{f: r[f] for f in FP}, label=label, confidence=conf,
                         defaulted=dflt))
    long_rows.append(dict(study_id=subj, reader_id=reader, session=session, detector_version=ver,
                          cand_uid=cand_uid(subj, r), list=lst, obj_id=oid, sheet_no=n, code=code,
                          label=label, coarse_call=call, confidence=conf, defaulted=dflt, blinded=1,
                          source="expert_sheet", timestamp=now))
out = f"{base}/review/{subj}_expert_{reader}.csv"
pd.DataFrame(rows, columns=["cand_id", "call"] + FP + ["label", "confidence", "defaulted"]).to_csv(out, index=False)
append_rows(f"{base}/labels/candidate_labels.csv", long_rows,
            ["study_id", "reader_id", "session", "cand_uid", "list"])
summ = pd.Series([r["call"] for r in rows]).value_counts().to_dict()
fine = pd.Series([r["label"] for r in long_rows]).value_counts().to_dict()
nx = sum(1 for r in long_rows if r["list"] == "excluded")
print(f"{subj}: expert {reader} -> {out}   {summ}")
print(f"{subj}: {len(long_rows)} labels ({nx} on excluded objects, "
      f"{sum(r['defaulted'] for r in long_rows)} defaulted to N) -> {base}/labels/candidate_labels.csv  {fine}")
if nx and any(r["coarse_call"] == "cSS" for r in long_rows if r["list"] == "excluded"):
    print(f"WARNING: the expert called {sum(1 for r in long_rows if r['list'] == 'excluded' and r['coarse_call'] == 'cSS')}"
          f" EXCLUDED object(s) cSS - an exclusion rule removed real cSS (see labels/candidate_labels.csv)")
if "--score" in sys.argv:
    css = sorted(r["cand_id"] for r in rows if r["call"] == "cSS")
    ich = sorted(r["cand_id"] for r in rows if r["call"] == "Near ICH")
    cmd = [sys.executable, f"{base}/scripts/mark_css.py", subj, ",".join(map(str, css)) or "none"]
    if ich: cmd += ["--ich", ",".join(map(str, ich))]
    subprocess.run(cmd, check=True)
```

## scripts/expert_sheet.py

```python
#!/usr/bin/env python3
"""Blinded picture sheet of all candidates for an expert reader -> review/ID_expert_sheet.pdf
usage: expert_sheet.py SUBJECT [--top N] [--seed 1] [--include-excluded] [--vessel-overlay]
  v4.16: every row also shows CORONAL and SAGITTAL views through the candidate (a vein that is a line
      on one plane is a dot on another; cSS stays a band along the cortex in all three).
  --vessel-overlay: outline the detector's vessel map (work/ID_vessels.nii.gz) in cyan on the minIP.
      Off by default: it is a detector output, so it can bias a blinded reader.
  --include-excluded (v4.15): the automatically EXCLUDED objects (review/ID_excluded.csv) are mixed in,
      unmarked, so the expert's calls show whether an exclusion rule removed real cSS and give
      labelled hard negatives (surface veins, artifacts, ...). Recommended for validation cases.

Each candidate gets a SHEET NUMBER in random order (the detector's rank is hidden, so the
order cannot bias the expert) and one row of pictures: SWI zoom (40 mm), the same with the
candidate outlined in red, and an 8 mm minIP slab (veins = continuous tubes). No features or
scores are printed. Page 1 explains the task; the last page is an answer grid.
The key (sheet number -> cand_id + fingerprint) is saved separately in
review/ID_expert_key.csv - do not give it to the expert.
Import the answers with expert_import.py."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, numpy as np, nibabel as nib, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, FP, LABEL_HELP

subj = sys.argv[1]
top = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else None
seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 1
base = _base()
swi = nib.as_closest_canonical(nib.load(f"{base}/data/{subj}_swi.nii"))
cnd = nib.as_closest_canonical(nib.load(f"{base}/review/{subj}_candidates.nii.gz"))
I = swi.get_fdata().astype(np.float32); C = cnd.get_fdata().astype(int)
vox = swi.header.get_zooms()[:3]
br = I > 0
lo, hi = np.percentile(I[br], [1, 99]) if br.any() else (I.min(), I.max())
df = pd.read_csv(f"{base}/review/{subj}_candidates.csv")
objs = [("candidate", int(c)) for c in df.cand_id.astype(int).tolist()[:top]]
X = None
if "--include-excluded" in sys.argv and os.path.exists(f"{base}/review/{subj}_excluded.csv"):
    xdf = pd.read_csv(f"{base}/review/{subj}_excluded.csv")
    if len(xdf):
        X = nib.as_closest_canonical(nib.load(f"{base}/review/{subj}_excluded.nii.gz")).get_fdata().astype(int)
        objs += [("excluded", int(e)) for e in xdf.excl_id]
order = np.random.default_rng(seed).permutation(len(objs))
sheet = {n + 1: objs[k] for n, k in enumerate(order)}          # sheet number -> (list, id)
rows = []
for n in sorted(sheet):
    lst, i = sheet[n]
    r = (df.set_index("cand_id") if lst == "candidate" else xdf.set_index("excl_id")).loc[i]
    rows.append(dict(sheet_no=n, cand_id=i if lst == "candidate" else None, **{f: r[f] for f in FP},
                     list=lst, obj_id=i))
key = pd.DataFrame(rows, columns=["sheet_no", "cand_id"] + FP + ["list", "obj_id"])
key.to_csv(f"{base}/review/{subj}_expert_key.csv", index=False)
mask_of = lambda n: (C == sheet[n][1]) if sheet[n][0] == "candidate" else (X == sheet[n][1])

HW = int(round(20 / vox[0]))                                   # 40 mm window
HWV = [int(round(20 / v)) for v in vox]                        # 40 mm in each axis (orthogonal views)
VS = None
if "--vessel-overlay" in sys.argv and os.path.exists(f"{base}/work/{subj}_vessels.nii.gz"):
    VS = nib.as_closest_canonical(nib.load(f"{base}/work/{subj}_vessels.nii.gz")).get_fdata() > 0
SLAB = max(1, int(round(4 / vox[2])))
PER_PAGE = 4
out = f"{base}/review/{subj}_expert_sheet.pdf"
with PdfPages(out) as pdf:
    f = plt.figure(figsize=(8.27, 11.69)); f.text(0.08, 0.92, f"cSS expert reading - {subj}", fontsize=16, weight="bold")
    f.text(0.06, 0.86, (
        "For each numbered candidate write one code (optionally ':1'-':5' = confidence,\n"
        "5 = certain; e.g. 12VS:4):\n\n"
        + "".join(f"   {c:<22s}{t}\n" for c, t in LABEL_HELP) + "\n"
        "Single letters C / V / N / H / U are enough if you prefer.\n\n"
        "Definition used (Charidimou et al., Neurology 2017): well-defined, homogeneous\n"
        "hypointense curvilinear signal on SWI outlining the outer surface of the cortex,\n"
        "within the adjacent subarachnoid space, or both.\n\n"
        "Each row: axial SWI (40 mm) | same, candidate outlined in red | 8 mm axial minIP\n"
        "(veins: continuous branching tubes) | coronal | sagittal (both outlined).\n"
        "A vein that is a line on one plane is usually a dot on another.\n"
        "Orientation: axial/coronal neurological (patient left on image left);\n"
        "sagittal: anterior on the right.\n\n"
        "Candidates are in RANDOM order; no scores are shown.\n"
        "Write the codes on the answer grid (last page)."), fontsize=9.5, va="top", family="monospace")
    f.text(0.06, 0.18, "Also: please mark any cSS you see that has NO number (outline it in freeview,\n"
           "or note slice and location) - it is needed to measure what the detector misses.", fontsize=10)
    pdf.savefig(f); plt.close(f)
    nums = sorted(sheet)
    for p0 in range(0, len(nums), PER_PAGE):
        f, axs = plt.subplots(PER_PAGE, 5, figsize=(8.27, 11.69))
        f.subplots_adjust(left=0.07, right=0.99, top=0.96, bottom=0.03, wspace=0.05, hspace=0.12)
        for r in range(PER_PAGE):
            for a in axs[r]: a.axis("off")
            if p0 + r >= len(nums): continue
            n = nums[p0 + r]; m = mask_of(n)
            if not m.any(): continue
            pts = np.argwhere(m); ci, cj = pts[:, 0].mean(), pts[:, 1].mean()
            kc = int(np.clip(np.bincount(pts[:, 2]).argmax(), 0, I.shape[2] - 1))
            i0, i1 = int(max(ci - HW, 0)), int(min(ci + HW, I.shape[0]))
            j0, j1 = int(max(cj - HW, 0)), int(min(cj + HW, I.shape[1]))
            sl = I[i0:i1, j0:j1, kc].T; ms = m[i0:i1, j0:j1, kc].T
            k0, k1 = max(kc - SLAB, 0), min(kc + SLAB + 1, I.shape[2])
            mip = np.where(I[i0:i1, j0:j1, k0:k1] > 0, I[i0:i1, j0:j1, k0:k1], hi).min(axis=2).T
            for c, (img, outline) in enumerate(((sl, False), (sl, True), (mip, True))):
                a = axs[r][c]; a.imshow(img, cmap="gray", origin="lower", vmin=lo, vmax=hi)
                if outline and ms.any():
                    src = m[i0:i1, j0:j1, k0:k1].any(axis=2).T if c == 2 else ms
                    a.contour(src.astype(float), levels=[0.5], colors="red", linewidths=0.9)
                if c == 2 and VS is not None and VS[i0:i1, j0:j1, k0:k1].any():
                    a.contour(VS[i0:i1, j0:j1, k0:k1].any(axis=2).T.astype(float), levels=[0.5],
                              colors="cyan", linewidths=0.5)
            # coronal (x-z) and sagittal (y-z) through the candidate centre, 40 mm windows, true aspect
            cjj, cii = int(round(cj)), int(round(ci))
            if not m[:, cjj, :].any(): cjj = int(np.bincount(pts[:, 1]).argmax())
            if not m[cii, :, :].any(): cii = int(np.bincount(pts[:, 0]).argmax())
            ka, kb = max(kc - HWV[2], 0), min(kc + HWV[2], I.shape[2])
            ia, ib = max(int(ci) - HWV[0], 0), min(int(ci) + HWV[0], I.shape[0])
            ja, jb = max(int(cj) - HWV[1], 0), min(int(cj) + HWV[1], I.shape[1])
            for c, (img, msk, asp) in ((3, (I[ia:ib, cjj, ka:kb], m[ia:ib, cjj, ka:kb], vox[2] / vox[0])),
                                       (4, (I[cii, ja:jb, ka:kb], m[cii, ja:jb, ka:kb], vox[2] / vox[1]))):
                a = axs[r][c]; a.imshow(img.T, cmap="gray", origin="lower", vmin=lo, vmax=hi, aspect=asp)
                if msk.any(): a.contour(msk.T.astype(float), levels=[0.5], colors="red", linewidths=0.9)
            axs[r][0].text(-0.06, 0.5, f"{n}", transform=axs[r][0].transAxes, fontsize=16, weight="bold",
                           ha="right", va="center")
            if r == 0 or p0 + r == p0:
                for c, t in enumerate(("axial", "outlined", "minIP 8 mm", "coronal", "sagittal")):
                    axs[r][c].set_title(t, fontsize=9)
        pdf.savefig(f); plt.close(f)
    # answer grid
    f = plt.figure(figsize=(8.27, 11.69)); f.text(0.08, 0.95, f"Answers - {subj}   (codes on page 1, e.g. C, VS:4)", fontsize=14, weight="bold")
    f.text(0.08, 0.925, "Reader: ____________________    Date: ____________", fontsize=10)
    cols = 4; rows = int(np.ceil(len(nums) / cols))
    for n in nums:
        c, r = (n - 1) // rows, (n - 1) % rows
        f.text(0.08 + c * 0.22, 0.89 - r * (0.85 / max(rows, 1)), f"{n:>3}  ______", fontsize=10, family="monospace")
    pdf.savefig(f); plt.close(f)
print(f"{subj}: expert sheet with {len(nums)} objects ({sum(1 for v in sheet.values() if v[0] == 'excluded')} "
      f"excluded mixed in) -> {out}\n"
      f"      key (do NOT give to the expert): {base}/review/{subj}_expert_key.csv")
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

## scripts/feature_report.py

```python
#!/usr/bin/env python3
"""Compare candidate features between calls (cSS vs Vein/Normal/Artifact) across all reviewed subjects.
usage: feature_report.py [--source expert|reader]
Default: for each subject the EXPERT reference (review/ID_expert_*.csv from expert_import.py) is used
when it exists, otherwise the reader's own calls (review/ID_calls.csv from review_css.py)."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, load_calls, match_calls
base = _base()
src_pref = sys.argv[sys.argv.index("--source") + 1] if "--source" in sys.argv else "expert"
subjects = sorted({os.path.basename(f).split("_calls.csv")[0] for f in glob.glob(f"{base}/review/*_calls.csv")} |
                  {os.path.basename(f).split("_expert_")[0] for f in glob.glob(f"{base}/review/*_expert_*.csv")
                   if not f.endswith("_expert_key.csv")})
rows = []; used = {}
for s in subjects:
    cp = f"{base}/review/{s}_candidates.csv"
    if not os.path.exists(cp): continue
    d = pd.read_csv(cp); calls = {}
    ex = sorted(f for f in glob.glob(f"{base}/review/{s}_expert_*.csv") if not f.endswith("_expert_key.csv"))
    if ex and src_pref == "expert":
        calls, dropped = match_calls(d, pd.read_csv(ex[0])); calls = calls or {}
        used[s] = "expert " + os.path.basename(ex[0]).split("_expert_")[1][:-4]
        if dropped: print(f"WARNING: {s}: {dropped} expert calls match no current candidate")
    else:
        calls, warn = load_calls(base, s, d)          # matched by fingerprint, not by rank
        if warn: print("WARNING:", warn)
        used[s] = "reader"
    if not calls: continue
    m = d.merge(pd.DataFrame(list(calls.items()), columns=["cand_id", "call"]), on="cand_id")
    m["subject"] = s; rows.append(m)
print("labels: " + ", ".join(f"{s}={u}" for s, u in used.items()))
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
         "mirror_dark_frac", "flair_csf_z", "parenchyma_frac", "css_evidence", "vein_evidence",
         "score_v3", "score_v4", "score"]
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

## scripts/label_missed.py

```python
#!/usr/bin/env python3
"""Record cSS that the detector did NOT show (needed to measure sensitivity, not just precision).
usage: label_missed.py SUBJECT READER [--mask path.nii.gz] [--session N]
The reader outlines every cSS focus that has no number on the expert sheet, in freeview on the SWI,
and saves it as labels/ID_missed_READER.nii.gz (any non-zero value; one connected blob per focus):
  freeview -v data/ID_swi.nii review/ID_candidates.nii.gz:colormap=lut  -> File > New Volume (template:
  the SWI) -> draw -> save as labels/ID_missed_READER.nii.gz
For each blob this script records where it went in the pipeline:
  shown       a current candidate covers >= 30 % of it (it was not missed; check the sheet call)
  excluded    an excluded object covers it -> the rule id and reason (review/ID_excluded.csv)
  dropped     removed by a generation gate -> the reason (work/ID_dropped.csv, v4.14)
  not generated  never a dark ridge component; v4.17 says why, from work/ID_gen.nii.gz (detect_css.py):
                 outside the search zone (skull-strip rim / white matter or deep-nucleus label / beyond the
                 pial band), not dark enough (z >= -2), dark but no ridge response, or no strong seed
-> labels/missed_lesions.csv (one row per blob; earlier rows of the same study/reader/session replaced)"""
import sys, os, datetime, numpy as np, pandas as pd, nibabel as nib
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, append_rows

if len(sys.argv) < 3: sys.exit(__doc__)
subj, reader = sys.argv[1], sys.argv[2]
base = _base()
arg = lambda f: sys.argv[sys.argv.index(f) + 1] if f in sys.argv else ""
session = int(arg("--session") or 1)
mp = os.path.expanduser(arg("--mask") or f"{base}/labels/{subj}_missed_{reader}.nii.gz")
if not os.path.exists(mp): sys.exit(f"no mask {mp} - draw it first (see usage)")
swi = nib.load(f"{base}/data/{subj}_swi.nii"); vox = np.array(swi.header.get_zooms()[:3], float)
M = nib.load(mp)
if M.shape[:3] != swi.shape[:3] or not np.allclose(M.affine, swi.affine, atol=1e-3):
    sys.exit("mask is not on the SWI grid - draw it on data/ID_swi.nii")
M = np.asarray(M.dataobj).reshape(swi.shape[:3]) > 0
ld = lambda p: nib.load(p).get_fdata().astype(int) if os.path.exists(p) else np.zeros(swi.shape[:3], int)
cand = ld(f"{base}/review/{subj}_candidates.nii.gz"); excl = ld(f"{base}/review/{subj}_excluded.nii.gz")
drop = ld(f"{base}/work/{subj}_dropped.nii.gz")
gen = ld(f"{base}/work/{subj}_gen.nii.gz")        # v4.17 generation-stage bits (all 0 for older detector runs)

def why_not_generated(m):
    """first generation stage the outline failed (bits from detect_css.py)"""
    g = gen[m]
    if not g.any(): return "unknown (no work/ID_gen.nii.gz - re-run detect_css.py v4.17)"
    if (g & 1).mean() < 0.5:
        rim_, wm_ = (g & 2 > 0).mean(), (g & 4 > 0).mean()
        if max(rim_, wm_) < 0.25: return "outside the search zone: beyond the pial band (-3.5..+4 mm)"
        return "outside the search zone: " + ("skull-strip rim (3 mm)" if rim_ >= wm_ else "white matter / deep-nucleus label")
    z_ = g[(g & 1) > 0]
    if not (z_ & 8).any(): return "not dark enough (z >= -2 everywhere)"
    if not ((z_ & 8 > 0) & (z_ & 32 > 0)).any(): return "dark but no ridge response (not line-shaped at the filter scales)"
    if not (g & 128).any(): return "no strong seed (needs z < z_thr AND a strong ridge for hysteresis)"
    return "generated, but no object covers it (split or merged away)"
xt = pd.read_csv(f"{base}/review/{subj}_excluded.csv").set_index("excl_id") \
    if os.path.exists(f"{base}/review/{subj}_excluded.csv") else pd.DataFrame()
dt = pd.read_csv(f"{base}/work/{subj}_dropped.csv").set_index("drop_id") \
    if os.path.exists(f"{base}/work/{subj}_dropped.csv") else pd.DataFrame()
seg_p = f"{base}/work/{subj}_seg_swispace.nii.gz"
seg = ld(seg_p)
lab, n = ndi.label(M, structure=np.ones((3, 3, 3), bool))
now = datetime.datetime.now().isoformat(timespec="seconds"); rows = []
top = lambda a, m: int(np.bincount(a[m][a[m] > 0]).argmax())
for i, sl in enumerate(ndi.find_objects(lab), 1):
    m = lab == i; near = ndi.binary_dilation(m, iterations=1)       # 1-voxel tolerance for blooming
    p = np.argwhere(m); c0 = p.mean(0)
    ext = float(np.sqrt((((p.max(0) - p.min(0) + 1) * vox) ** 2).sum()))
    cov = float((ndi.binary_dilation(cand > 0, iterations=1)[m]).mean())
    r = dict(study_id=subj, reader_id=reader, session=session, lesion_no=i,
             vox_i=int(round(c0[0])), vox_j=int(round(c0[1])), vox_k=int(round(c0[2])),
             volume_mm3=round(float(m.sum() * vox.prod()), 1), extent_mm=round(ext, 1),
             hemi="", region_label=0, fate="", obj_id=None, rule_version="", reason="",
             candidate_coverage=round(cov, 2), mask=os.path.basename(mp), timestamp=now)
    s_ = seg[near][seg[near] >= 1000]
    if s_.size:
        r["region_label"] = int(np.bincount(s_).argmax()); r["hemi"] = "L" if r["region_label"] < 2000 else "R"
    if cov >= 0.3:
        r["fate"], r["obj_id"] = "shown", top(cand, near)
    elif (excl[near] > 0).any():
        e = top(excl, near); r["fate"], r["obj_id"] = "excluded", e
        if e in xt.index:
            r["reason"] = xt.loc[e, "reason"]
            r["rule_version"] = xt.loc[e, "rule_version"] if "rule_version" in xt.columns else ""
    elif (drop[near] > 0).any():
        d = top(drop, near); r["fate"], r["obj_id"] = "dropped", d
        r["reason"] = dt.loc[d, "reason"] if d in dt.index else ""
    else:
        r["fate"] = "not generated"; r["reason"] = why_not_generated(m)
    rows.append(r)
# replace this study/reader/session's earlier rows (a redrawn mask may have fewer blobs)
path = f"{base}/labels/missed_lesions.csv"
if os.path.exists(path):
    old = pd.read_csv(path)
    old = old[~((old.study_id.astype(str) == subj) & (old.reader_id.astype(str) == reader) & (old.session == session))]
    old.to_csv(path, index=False)
append_rows(path, rows, ["study_id", "reader_id", "session", "lesion_no"])
f = pd.Series([r["fate"] for r in rows]).value_counts().to_dict()
print(f"{subj}: {n} missed-lesion outline(s) from {reader} -> {path}   {f}")
for r in rows:
    why = " ".join(x for x in (r["rule_version"], r["reason"]) if x)
    print(f"  #{r['lesion_no']}: {r['fate']}" + (f" ({why})" if why else "")
          + f"  extent {r['extent_mm']} mm  vox {r['vox_i']} {r['vox_j']} {r['vox_k']}")
if f.get("shown"):
    print("  note: 'shown' outlines are covered by a numbered candidate - check that number's call instead")
```

## scripts/make_synthetic.py

```python
#!/usr/bin/env python3
"""Insert synthetic cSS (and optional synthetic vein decoys) into a presumed cSS-negative host -> {HOST}S.
usage: make_synthetic.py HOST [n_lesions=8] [depth=0.6] [seed=1] [--tram 0.4] [--veins 4] [--radial 2]
                         [--surface-veins 2] [--legacy]

v3 (realistic lesions; after the P006 stress test showed v2 lesions were too big and blob-like, so
volume / branching / elongation separated them for the wrong reason):
  - thin CURVILINEAR lesions: a 1-voxel path traced along the cortex-CSF boundary (no flood fill),
    optionally thickened by 1 voxel (blooming), 5-30 mm long (log-uniform), on 2-5 consecutive slices
  - tram-track: with probability --tram the opposite bank of the same sulcus gets a matching line
  - darkness per lesion around `depth`, edge blur, and noise inside the lesion
  - --veins K: K tubular decoys (radius 0.5-0.9 mm, 15-40 mm) in the sulcal CSF, as dark as the
    lesions; written to work/HOSTS_veins.nii.gz (Cao et al. 2026: vessel decoys reduce false positives)
  - --radial K (v4.11): K transcortical / medullary vein decoys: straight tubes starting 2 mm out in the
    CSF and running perpendicular to the cortex 10-18 mm into the white matter (same veins file)
  - --surface-veins K (v4.16): the hardest real mimic - K veins running ALONG the cortical surface in the
    CSF within 1.5 mm of the cortex, 25-45 mm over many slices, 60 % with a side branch (8-15 mm) leaving
    10-20 mm along; same veins file; kinds in work/HOSTS_vein_kinds.csv (sulcal / radial / surface)
--physics (v4.18, cf. CenSynCMB 2026): instead of darkening voxels, every object gets a susceptibility (cSS
    hemosiderin +depth ppm, veins +0.45 ppm, calcium -1.0 ppm), the field is computed with the dipole kernel
    (css_common.susceptibility_signal, B0 along the scanner z axis), and the image is changed by intravoxel
    dephasing (blooming that grows with --te-ms and --b0 and depends on orientation), short-T2* signal loss inside
    the object, and the SWI phase mask^4 on the object's paramagnetic phase. Writes work/HOSTS_dphase.nii.gz.
    --b0 3.0 --te-ms 20 (defaults). An approximation for stress tests, not an SWI simulator.
--calcium K (with --physics): K calcification-like dots (r 0.8-1.5 mm, diamagnetic) at the cortex - a cSS /
    microbleed mimic that SWI magnitude alone cannot tell apart (phase can); kind "calcium" in the veins file.
--legacy reproduces the v2 generator (flood-filled, thicker lesions) for comparison."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib
from scipy import ndimage as ndi

flag = lambda f, d: type(d)(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d
TRAM, NVEIN, NRAD, LEGACY = flag("--tram", 0.4), flag("--veins", 0), flag("--radial", 0), "--legacy" in sys.argv
NSURF = flag("--surface-veins", 0)
PHYS, B0, TE_MS, NCA = "--physics" in sys.argv, flag("--b0", 3.0), flag("--te-ms", 20.0), flag("--calcium", 0)
# positional values given after a flag (e.g. "--veins 4") must not be read as lesion arguments
pos = [a for i, a in enumerate(sys.argv[1:], 1) if not a.startswith("--") and not sys.argv[i - 1] in ("--tram", "--veins", "--radial", "--surface-veins", "--b0", "--te-ms", "--calcium")]
host = pos[0]; N = int(pos[1]) if len(pos) > 1 else 8
depth = float(pos[2]) if len(pos) > 2 else 0.6; seed = int(pos[3]) if len(pos) > 3 else 1
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
pts_all = pts_all[(pts_all[:, 2] >= 3) & (pts_all[:, 2] < nz - 6)]
if len(pts_all) == 0: sys.exit("no admissible cortical boundary found")
NB8 = [(dx, dy) for dx in (-1, 0, 1) for dy in (-1, 0, 1) if dx or dy]

def trace(b, start, L_mm, d0=None):
    """1-voxel path along the boundary b (2-D), preferring to keep its direction (no branches)"""
    path, seen, cur, d, length = [start], {start}, start, d0, 0.0
    while length < L_mm:
        nxt = [(cur[0] + dx, cur[1] + dy) for dx, dy in NB8
               if 0 <= cur[0] + dx < b.shape[0] and 0 <= cur[1] + dy < b.shape[1]
               and b[cur[0] + dx, cur[1] + dy] and (cur[0] + dx, cur[1] + dy) not in seen]
        if not nxt: break
        if d is None: n = nxt[rng.integers(len(nxt))]
        else: n = max(nxt, key=lambda q: (q[0] - cur[0]) * d[0] + (q[1] - cur[1]) * d[1])
        d = (n[0] - cur[0], n[1] - cur[1])
        length += float(np.hypot(d[0] * vox[0], d[1] * vox[1]))
        path.append(n); seen.add(n); cur = n
    return path, length, d

def opposite_bank(path, k):
    """boundary voxels across the sulcal CSF from the path (1.5-6 mm away, CSF in between)"""
    b = boundary[:, :, k]; P = np.array(path)
    cand = np.argwhere(b)
    if not len(cand): return None
    dd = np.sqrt((((cand[:, None, :] - P[None, :, :]) * vox[:2]) ** 2).sum(-1)).min(1)
    cand = cand[(dd > 1.5) & (dd < 6.0)]
    ok = [q for q in map(tuple, cand)
          if sulcal[(q[0] + P[np.argmin(((P - q) ** 2).sum(1))][0]) // 2,
                    (q[1] + P[np.argmin(((P - q) ** 2).sum(1))][1]) // 2, k]]
    return ok[rng.integers(len(ok))] if ok else None

truth, seeds, tries = np.zeros(I.shape, np.int16), [], 0
ldepth = {}
dark = np.zeros(I.shape, np.float32)
while len(seeds) < N and tries < 8000:
    tries += 1
    p = pts_all[rng.integers(len(pts_all))]
    if any(np.linalg.norm((p - q) * vox) < 25 for q in seeds): continue
    lesion, ok = np.zeros(I.shape, bool), True
    if LEGACY:                                   # v2: flood fill, then 1-2 dilations
        L_mm, nsl, cur = rng.uniform(10, 30), int(rng.integers(2, 5)), p[:2]
        for dk in range(nsl):
            k = p[2] + dk; b = boundary[:, :, k]; cand = np.argwhere(b)
            if k >= nz or len(cand) == 0: ok = False; break
            d = np.linalg.norm((cand - cur) * vox[:2], axis=1); j = d.argmin()
            if d[j] > 4: ok = False; break
            start = tuple(cand[j]); cur = cand[j]
            target, seen, queue = int(L_mm / vox[0]), {start}, [start]
            while queue and len(seen) < target:
                x, y = queue.pop(0)
                for dx, dy in NB8:
                    q = (x + dx, y + dy)
                    if 0 <= q[0] < b.shape[0] and 0 <= q[1] < b.shape[1] and b[q] and q not in seen:
                        seen.add(q); queue.append(q)
            if len(seen) < 0.6 * target: ok = False; break
            for q in seen: lesion[q[0], q[1], k] = True
        if not ok: continue
        lesion = ndi.binary_dilation(lesion, structure=s2, iterations=int(rng.integers(1, 3))) & brain
    else:                                        # v3: thin traced line(s)
        L_mm = float(np.exp(rng.uniform(np.log(5), np.log(30))))
        nsl = int(rng.integers(2, 6)); cur, d0 = tuple(p[:2]), None
        tram = rng.random() < TRAM; opp = None
        for dk in range(nsl):
            k = p[2] + dk; b = boundary[:, :, k]; cand = np.argwhere(b)
            if k >= nz or not len(cand): ok = False; break
            dd = np.linalg.norm((cand - np.array(cur)) * vox[:2], axis=1); j = dd.argmin()
            if dd[j] > 3: ok = False; break
            path, length, d0 = trace(b, tuple(cand[j]), L_mm, d0 if dk else None)
            if length < 0.6 * L_mm: ok = False; break
            for q in path: lesion[q[0], q[1], k] = True
            cur = tuple(cand[j])
            if tram:
                o = opposite_bank(path, k) if opp is None else None
                start = o if o is not None else opp
                if start is not None:
                    bo = boundary[:, :, k]; cand_o = np.argwhere(bo)
                    jo = np.linalg.norm((cand_o - np.array(start)) * vox[:2], axis=1).argmin()
                    po, _, _ = trace(bo, tuple(cand_o[jo]), L_mm)
                    for q in po: lesion[q[0], q[1], k] = True
                    opp = tuple(cand_o[jo])
        if not ok: continue
        if rng.random() < 0.5:                   # blooming: 1-voxel in-plane thickening, stays at the surface
            lesion = ndi.binary_dilation(lesion, structure=s2) & (cortex | sulcal) & brain
    if (truth[lesion] > 0).any() or not lesion.any(): continue
    seeds.append(p); truth[lesion] = len(seeds)
    d_i = float(np.clip(depth + rng.uniform(-0.1, 0.1), 0.05, 0.9)); ldepth[len(seeds)] = d_i
    w = ndi.gaussian_filter(lesion.astype(np.float32), sigma=(rng.uniform(0.4, 0.8),) * 2 + (0,))
    dark = np.maximum(dark, np.clip(w * 2.0, 0, 1) * d_i)

# ---- vein decoys: tubes running in the sulcal CSF
veins = np.zeros(I.shape, np.int16)
if NVEIN and not LEGACY:
    mid = sulcal & (ndi.distance_transform_edt(sulcal, sampling=vox) >= 0.7) & (dist > 6.0)
    in_csf = ndi.binary_dilation(sulcal, structure=np.ones((3, 3, 3), bool))   # decoys stay in sulcal CSF
    vp = np.argwhere(mid); nv = 0; vt = 0
    while nv < NVEIN and vt < 2000 and len(vp):
        vt += 1
        q = vp[rng.integers(len(vp))].astype(float)
        if (truth[tuple(q.astype(int))] > 0) or any(np.linalg.norm((q - s) * vox) < 15 for s in seeds): continue
        dvec = rng.normal(size=3); dvec[2] *= rng.uniform(0.2, 1.0); dvec /= np.linalg.norm(dvec * vox) + 1e-9
        L, r, pts, pnt = rng.uniform(15, 40), rng.uniform(0.5, 0.9), [], q.copy()
        for _ in range(int(L / 0.5)):
            dvec += rng.normal(scale=0.05, size=3) / vox; dvec /= np.linalg.norm(dvec * vox) + 1e-9
            pnt = pnt + dvec * 0.5
            ii = np.round(pnt).astype(int)
            if (ii < 0).any() or (ii >= np.array(I.shape)).any() or not in_csf[tuple(ii)]: break
            pts.append(ii)
        if len(pts) * 0.5 < 10: continue
        m = np.zeros(I.shape, bool); m[tuple(np.array(pts).T)] = True
        m = (ndi.distance_transform_edt(~m, sampling=vox) <= r) & brain
        if (truth[m] > 0).any() or (veins[m] > 0).any(): continue
        nv += 1; veins[m] = nv
        dark = np.maximum(dark, m * float(np.clip(depth + rng.uniform(-0.1, 0.1), 0.05, 0.9)))
n_sulc = int(veins.max())

# ---- transcortical / medullary vein decoys: straight tubes perpendicular to the cortex into the WM
if NRAD and not LEGACY:
    wm = np.isin(seg, [2, 41])
    g = np.gradient(ndi.gaussian_filter(ndi.distance_transform_edt(~wm, sampling=vox), 1.0), *vox)
    nv, vt = int(veins.max()), 0; n0 = nv
    while nv - n0 < NRAD and vt < 2000:
        vt += 1
        q = pts_all[rng.integers(len(pts_all))].astype(float)
        if any(np.linalg.norm((q - s_) * vox) < 15 for s_ in seeds) or veins[tuple(q.astype(int))]: continue
        u = -np.array([gg[tuple(q.astype(int))] for gg in g]); nu = np.linalg.norm(u)
        if nu < 0.3: continue
        u = u / nu                                         # mm direction toward the white matter
        L, r = rng.uniform(10, 18), rng.uniform(0.5, 0.9)
        tt = np.arange(-2.0, 3.0 + L, 0.3)                 # 2 mm outside, ~3 mm cortex, L mm into WM
        ii = np.round(q[None, :] + tt[:, None] * u[None, :] / vox).astype(int)
        if (ii < 0).any() or (ii >= np.array(I.shape)).any(): continue
        if wm[tuple(ii[tt > 4.0].T)].mean() < 0.6 or not brain[tuple(ii.T)].all(): continue
        m = np.zeros(I.shape, bool); m[tuple(ii.T)] = True
        m = (ndi.distance_transform_edt(~m, sampling=vox) <= r) & brain
        if (truth[m] > 0).any() or (veins[m] > 0).any(): continue
        nv += 1; veins[m] = nv
        dark = np.maximum(dark, m * float(np.clip(depth + rng.uniform(-0.1, 0.1), 0.05, 0.9)))

# ---- surface-following vein decoys (v4.16): direction-persistent walk that stays in the CSF next to the cortex
kinds = {v: ("sulcal" if v <= n_sulc else "radial") for v in range(1, int(veins.max()) + 1)}
if NSURF and not LEGACY:
    near_ctx = brain & ~cortex & ~np.isin(seg, [2, 41]) & \
        (ndi.distance_transform_edt(~cortex, sampling=vox) <= 1.5) & (dist > 6.0)
    sp = np.argwhere(near_ctx); n0 = int(veins.max()); ns = vt = 0
    def walk(start, d, L):
        pts, pnt = [], start.astype(float)
        for _ in range(int(L / 0.4)):
            for t in range(25):                      # keep direction, else turn to stay next to the cortex
                dv = d + rng.normal(scale=0.08 if t == 0 else 0.5, size=3) / vox
                dv /= np.linalg.norm(dv * vox) + 1e-9
                q = pnt + dv * 0.4; ii = np.round(q).astype(int)
                if (ii >= 0).all() and (ii < np.array(I.shape)).all() and near_ctx[tuple(ii)]:
                    break
            else:
                break
            d, pnt = dv, q; pts.append(np.round(pnt).astype(int))
        return pts, d
    while ns < NSURF and vt < 3000 and len(sp):
        vt += 1
        q = sp[rng.integers(len(sp))]
        if truth[tuple(q)] or veins[tuple(q)] or any(np.linalg.norm((q - s_) * vox) < 15 for s_ in seeds): continue
        d0 = rng.normal(size=3); d0[2] *= rng.uniform(0.3, 1.0); d0 /= np.linalg.norm(d0 * vox) + 1e-9
        L = rng.uniform(25, 45); pts, _ = walk(q, d0, L)
        if len(pts) * 0.4 < 0.7 * L: continue
        if rng.random() < 0.6:                       # side branch 10-20 mm along the vein
            k = min(int(rng.uniform(10, 20) / 0.4), len(pts) - 2)
            dd = np.array(pts[k + 1], float) - np.array(pts[k], float)
            br = np.cross(dd, rng.normal(size=3)); br /= np.linalg.norm(br * vox) + 1e-9
            bp, _ = walk(np.array(pts[k]), 0.5 * dd / (np.linalg.norm(dd * vox) + 1e-9) + br, rng.uniform(8, 15))
            pts += bp
        m = np.zeros(I.shape, bool); m[tuple(np.array(pts).T)] = True
        m = (ndi.distance_transform_edt(~m, sampling=vox) <= rng.uniform(0.5, 0.8)) & brain
        if (truth[m] > 0).any() or (veins[m] > 0).any(): continue
        ns += 1; veins[m] = n0 + ns; kinds[n0 + ns] = "surface"
        dark = np.maximum(dark, m * float(np.clip(depth + rng.uniform(-0.1, 0.1), 0.05, 0.9)))

# ---- calcification-like dots (v4.18, --physics only): diamagnetic, at the cortex
if NCA and PHYS and not LEGACY:
    n0, nc, ct = int(veins.max()), 0, 0
    while nc < NCA and ct < 2000:
        ct += 1
        q = pts_all[rng.integers(len(pts_all))]
        if any(np.linalg.norm((q - s_) * vox) < 12 for s_ in seeds) or veins[tuple(q)]: continue
        rr = rng.uniform(0.8, 1.5)
        m = np.zeros(I.shape, bool); m[tuple(q)] = True
        m = (ndi.distance_transform_edt(~m, sampling=vox) <= rr) & brain
        if (truth[m] > 0).any() or (veins[m] > 0).any(): continue
        nc += 1; veins[m] = n0 + nc; kinds[n0 + nc] = "calcium"

if PHYS:
    from css_common import susceptibility_signal
    te = TE_MS / 1000.0
    css_fill = np.zeros(I.shape, np.float32); chi = np.zeros(I.shape, np.float32)
    for i, d_i in ldepth.items():                       # thin hemosiderin layer, partial volume -> blurred fill
        w = np.clip(2.0 * ndi.gaussian_filter((truth == i).astype(np.float32), sigma=(0.5, 0.5, 0)), 0, 1)
        css_fill = np.maximum(css_fill, w); chi += w * d_i * 1.0           # ppm: depth 0.6 -> 0.6 ppm
    is_ca = np.isin(veins, [v for v, k_ in kinds.items() if k_ == "calcium"])
    ves_fill = ((veins > 0) & ~is_ca).astype(np.float32)
    chi += 0.45 * ves_fill - 1.0 * is_ca
    ph, att = susceptibility_signal(chi, vox, img.affine, B0, te)
    # signal inside the object relative to tissue: exp(-TE * dR2*), dR2* 60 /s hemosiderin, 30 /s venous blood
    inside = 1 - css_fill * (1 - np.exp(-te * 60.0)) - ves_fill * (1 - np.exp(-te * 30.0)) - is_ca * 0.5
    swi_mask = np.where(ph > 0, np.clip(1 - ph / np.pi, 0, 1), 1.0) ** 4    # SWI phase mask (paramagnetic +)
    Is = I * att * np.clip(inside, 0, 1) * swi_mask
    dark = 1 - Is / np.maximum(I, 1e-6)
    nib.save(nib.Nifti1Image(ph, img.affine), f"{base}/work/{host}S_dphase.nii.gz")
else:
    Is = I * (1 - dark)
Is += (rng.normal(0, 0.35 * sd_c, I.shape) * (dark > 0.05)).astype(np.float32)
Is[brain] = np.maximum(Is[brain], 1.0)        # keep the brain mask unchanged
Is[~brain] = I[~brain]

out = nib.Nifti1Image(Is.astype(np.float32), img.affine, img.header); out.set_data_dtype(np.float32)
nib.save(out, f"{base}/data/{host}S_swi.nii")
nib.save(nib.Nifti1Image(truth, img.affine), f"{base}/work/{host}S_truth.nii.gz")
nib.save(nib.Nifti1Image(veins, img.affine), f"{base}/work/{host}S_veins.nii.gz")
import pandas as pd
pd.DataFrame(sorted(kinds.items()), columns=["vein_id", "kind"]).to_csv(f"{base}/work/{host}S_vein_kinds.csv", index=False)

lut = {}
lutp = os.path.join(os.environ.get("FREESURFER_HOME", ""), "FreeSurferColorLUT.txt")
if os.path.exists(lutp):
    for line in open(lutp):
        t = line.split()
        if len(t) > 2 and t[0].isdigit(): lut[int(t[0])] = t[1]
print(f"{host}S: inserted {len(seeds)} synthetic cSS lesions ({'legacy' if LEGACY else 'v3 thin'}, depth ~{depth}, "
      f"seed {seed}){f', {int(veins.max())} vein / mimic decoys' if veins.max() else ''}"
      + (f"; physics model B0 {B0:g} T, TE {TE_MS:g} ms" if PHYS else ""))
for i in range(1, len(seeds) + 1):
    m = truth == i; s = seg[m][seg[m] >= 1000]
    lab = int(np.bincount(s).argmax()) if s.size else 0
    ks = np.unique(np.argwhere(m)[:, 2]); c0 = np.argwhere(m).mean(0).round().astype(int)
    print(f"  lesion {i}: {lut.get(lab, lab)!s:30s} slices {ks.min()}-{ks.max()}  vox {c0[0]} {c0[1]} {c0[2]}"
          f"  {m.sum() * np.prod(vox):.0f} mm3")
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

## scripts/review_css.py

```python
#!/usr/bin/env python3
"""Interactive cSS review: shows each suspected candidate, you agree or disagree.
usage: review_css.py SUBJECT [--top N] [--redo] [--by-tier]
  --by-tier (v4.18, opt-in): tier-1 candidates first, then tier 2 (review_tier, recorded by detect_css.py;
  not validated). Every candidate is still shown; with --top N the cut applies to the reordered list.
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
Right column (only with run_css.sh --phase): filtered SWI phase - calcium has the opposite sign to
veins and blood products (AJNR 2016 mimic list).
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
_qcp = f"{base}/review/{subj}_qc.json"          # v4.14 case QC (detect_css.py)
if os.path.exists(_qcp):
    import json as _json
    _qc = _json.load(open(_qcp))
    if _qc.get("status") != "ok":
        print(f"WARNING: case QC {_qc.get('status')}: {'; '.join(_qc.get('reasons', []))} - "
              f"anatomy-based features and exclusions are less reliable for this scan")
if "--by-tier" in sys.argv and "review_tier" in df.columns:
    ids = df.sort_values(["review_tier", "cand_id"], kind="stable").cand_id.astype(int).tolist()[:top]
    print(f"order: tier 1 first ({int((df.review_tier == 1).sum())}), then tier 2 ({int((df.review_tier == 2).sum())})")
else:
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
# optional SWI phase (run_css.sh --phase): high-pass filtered phase separates paramagnetic blood
# products (same sign as veins) from diamagnetic CALCIUM (opposite sign) - an AJNR 2016 cSS mimic.
# Sign conventions differ by vendor (Vaccarino et al. 2025), so compare with a vein on the same image.
PH = None
_php = f"{base}/data/{subj}_phase.nii.gz"
if os.path.exists(_php) and os.path.getmtime(_php) < os.path.getmtime(f"{base}/data/{subj}_swi.nii"):
    print("WARNING: phase image is older than the SWI - not shown (re-import it with run_css.sh --phase)")
elif os.path.exists(_php):
    from nibabel.processing import resample_from_to
    _ph = nib.load(_php)
    if _ph.shape[:3] != swi.shape[:3] or not np.allclose(_ph.affine, swi.affine, atol=1e-3):
        _ph = resample_from_to(_ph, swi, order=0)
    else:
        _ph = nib.as_closest_canonical(_ph)
    PH = _ph.get_fdata().astype(np.float32)
    _lo, _hi = float(PH.min()), float(PH.max())
    if _hi - _lo > 2 * np.pi + 0.5:                 # scanner units (e.g. -4096..4095) -> radians
        PH = (PH - _lo) / (_hi - _lo) * 2 * np.pi - np.pi
from scipy import ndimage as _ndi
def hp_phase(sl2d):
    """homodyne high-pass of one phase slice: angle(z * conj(lowpass z)), wrap-safe"""
    z = np.exp(1j * sl2d); sg = 4.0 / vox[0]
    lp = _ndi.gaussian_filter(z.real, sg) + 1j * _ndi.gaussian_filter(z.imag, sg)
    return np.angle(z * np.conj(lp))
SLAB = max(1, int(round(4 / vox[2]))) # minIP slab half-thickness (~8 mm total)
def num(v):
    try: return float(v)
    except (TypeError, ValueError): return float("nan")
flag = lambda v: num(v) > 0          # NaN / missing column -> False

class Reviewer:
    def __init__(self):
        self.pos = next((n for n, c in enumerate(ids) if c not in calls), 0)
        self.shift = 0
        ncol = 5 if PH is not None else 4
        self.fig = plt.figure(figsize=(14 + 3.5 * (ncol - 4), 8.2))
        gs = self.fig.add_gridspec(2, ncol, left=0.02, right=0.98, top=0.85, bottom=0.14, wspace=0.05, hspace=0.12)
        self.ph = [self.fig.add_subplot(gs[r, 4]) for r in range(2)] if PH is not None else None
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
        if self.ph is not None:
            hp = hp_phase(PH[:, :, kc])[i0:i1, j0:j1]; br2 = I[i0:i1, j0:j1, kc] > 0
            v = float(np.percentile(np.abs(hp[br2]), 98)) if br2.any() else 1.0
            for row in (0, 1):
                a = self.ph[row]; a.clear()
                a.imshow(np.where(br2, hp, 0).T, cmap="gray", origin="lower", vmin=-v, vmax=v)
                if row == 1 and m[i0:i1, j0:j1, kc].any():
                    a.contour(m[i0:i1, j0:j1, kc].T.astype(float), levels=[0.5], colors="red", linewidths=1.2)
                a.set_xticks([]); a.set_yticks([])
            self.ph[0].set_title("filtered PHASE (centre slice)\ncalcium = opposite sign to veins", fontsize=9)
        flags = [n for n, f in (("artifact zone", "artifact_zone"), ("midline", "midline_zone"), ("vein-like", "vein_like"),
                                ("INFRATENTORIAL", "infratentorial"), ("FLAIR-bright CSF: acute cSAH?", "flair_bright"),
                                ("FLAIR-bright cortex: cortical vein thrombosis?", "flair_ctx_bright"),
                                ("microbleed-like: small, >=half in parenchyma", "cmb_like")) if flag(r.get(f, 0))]
        if flag(r.get("near_ich_suggest", 0)):
            flags.append(f"ICH {num(r.get('ich_dist_mm')):.0f} mm away - press i if contiguous")
        feat = ""
        if "pial_dist_mm" in r:
            shape = "plate" if num(r.tube_ratio) < 0.4 else ("tube" if num(r.tube_ratio) > 0.5 else "mixed")
            tram = num(r.tram_frac)
            ev = (f"   |   evidence: cSS {num(r.css_evidence):.2f}  vein {num(r.vein_evidence):.2f}"
                  if "css_evidence" in r else "")
            feat = (f"\npial {num(r.pial_dist_mm):+.1f} mm   shape {shape} ({num(r.tube_ratio):.2f})   "
                    f"follows surface {num(r.surface_alignment):.2f}   "
                    f"tram-track {'n/a' if np.isnan(tram) else f'{tram:.0%}'}   vein tree {num(r.vein_tree_mm):.0f} mm   "
                    f"mirror dark {num(r.mirror_dark_frac):.0%}" + ev)
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

## scripts/rule_test.py

```python
#!/usr/bin/env python3
"""Which single-feature exclusion rule would remove false positives WITHOUT losing known cSS?
usage: rule_test.py --pos P006 --old OLD_candidates.nii.gz --ids 3,4,5,13,18,19 --neg P011[,P010]
                    [--synthetic results/stress_TAG_features.csv] [--top 12]
  --pos/--old/--ids  known cSS: lesions of an older candidate map (as in check_known.py), matched to
                     the CURRENT candidates of that subject
  --neg              cSS-negative subjects: every current candidate there is a false positive
  --synthetic        optional stress-test feature table: rules are also checked on synthetic lesions
                     (is_lesion == 1); synthetic features are partly artifacts of the generator
For each feature and direction the threshold is put just beyond the most extreme known lesion (plus a
safety margin), so every known lesion is kept; the table shows how many false positives that removes.
With only a few known lesions this is a SHORTLIST for discussion, not a validated rule."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, numpy as np, pandas as pd, nibabel as nib
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base

base = _base()
arg = lambda f, d="": sys.argv[sys.argv.index(f) + 1] if f in sys.argv else d
if not (arg("--pos") and arg("--old") and arg("--ids") and arg("--neg")): sys.exit(__doc__)
pos_s, top = arg("--pos"), int(arg("--top", "12"))
FEATS = ["darkness_z", "volume_mm3", "elongation", "surface_contact", "cortex_frac", "cortex_dist_mm",
         "surface_gradient", "pial_dist_mm", "bank_frac", "tube_ratio", "surface_alignment", "tram_frac",
         "vein_tree_mm", "mirror_dark_frac", "edge_contrast", "parenchyma_frac", "extent_mm",
         "vessel_ext_mm", "vessel_wm_mm", "vessel_run_mm", "axis_normal", "css_evidence", "vein_evidence",
         "score_v4"]

# known lesions -> current candidates of the positive subject
old = nib.load(os.path.expanduser(arg("--old"))).get_fdata().astype(int)
cand = nib.load(f"{base}/review/{pos_s}_candidates.nii.gz").get_fdata().astype(int)
if old.shape != cand.shape: sys.exit("grids differ - re-import the same SWI")
dp = pd.read_csv(f"{base}/review/{pos_s}_candidates.csv")
known = {}
for i in [int(x) for x in arg("--ids").split(",") if x.strip()]:
    v = cand[ndi.binary_dilation(old == i, iterations=2)]; v = v[v > 0]
    if len(v): known[i] = int(np.bincount(v).argmax())          # the candidate covering it most
    else: print(f"known lesion #{i}: not a current candidate (excluded or not detected) - not used")
if not known: sys.exit("no known lesion is a current candidate")
P = dp[dp.cand_id.isin(set(known.values()))]
O = dp[~dp.cand_id.isin(set(known.values()))]          # other candidates of the positive subject
N = pd.concat([pd.read_csv(f"{base}/review/{s}_candidates.csv").assign(subject=s)
               for s in arg("--neg").split(",")], ignore_index=True)
S = None
if arg("--synthetic"):
    S = pd.read_csv(os.path.expanduser(arg("--synthetic"))); S = S[S.is_lesion == 1]
F = [f for f in FEATS if f in dp.columns and f in N.columns]

print(f"known cSS ({pos_s}): {len(P)} candidates {sorted(known.items())}  (old id: current id)")
print(f"false positives: {len(N)} candidates of {arg('--neg')}" +
      (f";  synthetic lesions: {len(S)}" if S is not None else "") +
      f";  other {pos_s} candidates (unlabelled): {len(O)}\n")
print("Feature values of the known lesions vs the false positives (median [min-max]):")
for f in F:
    p_, n_ = P[f].dropna(), N[f].dropna()
    if not len(p_) or not len(n_): continue
    print(f"  {f:18s} known {p_.median():7.2f} [{p_.min():7.2f} {p_.max():7.2f}]   "
          f"FP {n_.median():7.2f} [{n_.min():7.2f} {n_.max():7.2f}]")

rows = []
for f in F:
    p_ = P[f].dropna()
    if len(p_) < len(P) or not N[f].notna().any(): continue      # a rule must be defined for every known lesion
    allv = pd.concat([p_, N[f].dropna()]); margin = 0.1 * (allv.max() - allv.min())
    for op, t in ((">", p_.max() + margin), ("<", p_.min() - margin)):
        hit = (lambda x: x > t) if op == ">" else (lambda x: x < t)
        nrem = int(hit(N[f].dropna()).sum())
        if nrem == 0: continue
        rows.append(dict(rule=f"exclude if {f} {op} {t:.2f}", fp_removed=nrem,
                         fp_pct=round(100 * nrem / len(N)),
                         other_pos_removed=int(hit(O[f].dropna()).sum()),
                         synthetic_lost=(f"{int(hit(S[f].dropna()).sum())}/{len(S)}" if S is not None and f in S else "-"),
                         margin_to_known=round(margin, 2)))
print("\nRules that keep ALL known lesions (threshold = most extreme known lesion + 10 % of the range):")
if rows:
    R = pd.DataFrame(rows).sort_values("fp_removed", ascending=False).head(top)
    print(R.to_string(index=False))
else:
    print("  none - every feature of the false positives overlaps the known lesions")
print("\nRead with care: few known lesions -> a rule can look safe by chance. Prefer rules that also lose "
      "few synthetic lesions and that make sense radiologically; confirm on the next patients.")
```

## scripts/score_css.py

```python
#!/usr/bin/env python3
"""Score accepted candidates.  usage: score_css.py SUBJECT [--truth] [--explain]
--explain: print which sulci each focus lines and why sulci count as adjacent (shared gyrus / touch)
v2: adds van Harten-style seeded region growing - accepted candidates grow into connected
dark voxels (z<-1.5, inside the search zone, max 5 mm from the candidate) to give the
full lesion extent and a continuous cSS volume (no 0-4 ceiling effect).
Candidates marked near_ich=1 (siderosis connected to a lobar ICH) are excluded and reported.
v4: SULCAL scoring when work/ID_a2009s_swispace.nii.gz exists (recon-all + prep_anat.sh):
  each accepted focus is assigned to the Destrieux sulci it lines (within 4 mm; >=15 % of the focus
  or >=20 mm3 of it); two sulci are
  "immediately adjacent" if they touch or border the same gyrus. Per hemisphere (Charidimou):
  0 none; 1 = one sulcus or <=3 adjacent sulci; 2 = >=2 non-adjacent or >3 sulci.
  STRIVE-2 category: focal = 1-3 sulci, disseminated = >3 sulci.
  Without Destrieux labels the v3 Euclidean approximation is used (3 mm foci, 10 mm adjacency).
Infratentorial candidates (classical superficial siderosis pattern) are reported, not scored.
v4.13 (van Harten 2023 Discussion): cortical surface area covered by cSS (mm2) and % of each hemisphere's
cortical surface - less blooming- and protocol-dependent than volume.
v4.3 ICH rule (Charidimou et al., Neurology 2017;89:2128): cSS "contiguous or potentially anatomically
connected with any lobar ICH" is not scored; cSS must be separated from any lobar ICH by >=3
unaffected sulci, or by >=2 (at multiple axial levels) if the haematoma has no superficial path along
the convexity. With Destrieux labels + an ICH mask, unaffected sulci between each accepted focus and
the ICH are counted on the sulcal adjacency graph: <2 -> excluded (reader-drawn ICH mask) or warned
(automatic mask); exactly 2 -> kept but flagged for the reader to check the 2-sulci conditions.
v4.4 Boston criteria v2.0 (Charidimou et al., Lancet Neurol 2022;21:714) count cSS by GYRI: a single
focus (even extending to a second adjacent gyrus) = 1 haemorrhagic lesion; multifocal cSS (gyri
separated by uninvolved areas, or >=3 adjacent gyri) = >=2 lesions. Reported as boston2_css_lesions."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, json, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import load_lut, fresh

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
A2 = nib.load(a2_p).get_fdata().astype(np.int32) if fresh(a2_p, f"{base}/data/{subj}_swi.nii") else None
if os.path.exists(a2_p) and A2 is None:          # v4.17: Destrieux labels from before the SWI was re-imported
    print("WARNING: Destrieux labels are older than the SWI - ignored (re-run prep_anat.sh); Euclidean scoring")
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
    # a sulcus counts when it holds >=15 % of the focus OR >=20 mm3 of it (v4.8: merged foci that
    # run from one sulcus over the crown into the next must keep both sulci)
    keep = (cnt >= 3) & ((cnt >= 0.15 * len(codes)) | (cnt * vmm3 >= 20.0))
    return [int(c + lo) for c in np.nonzero(keep)[0]]

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
        fs = {i: cand_sulci.get(i) or sulci_of(masks[i], h) for i in ids}
        S = sorted(set(c for i in ids for c in fs[i]))
        if "--explain" in sys.argv and S:
            nm = lambda c: lut.get(c, str(c)).replace("ctx_lh_", "").replace("ctx_rh_", "")
            print(f"\n[{h}] foci -> sulci")
            for i in ids:
                print(f"   focus #{i}: {', '.join(nm(c) for c in fs[i]) or '-'}")
            print(f"[{h}] sulcus pairs: adjacency link | closest foci (mm)")
            for a_ in S:
                for b_ in S:
                    if a_ >= b_: continue
                    ta, tb = touching(a_), touching(b_)
                    link = "touch" if (b_ in ta or a_ in tb) else \
                        ", ".join(nm(g) for g in sorted(ta & tb) if is_gyr(g)) or "NOT adjacent"
                    dd = [D[(min(i, j), max(i, j))] for i in ids for j in ids
                          if i != j and a_ in fs[i] and b_ in fs[j]]
                    dmin = f"{min(dd):.0f}" if dd else ("same focus" if any(a_ in fs[i] and b_ in fs[i] for i in ids) else "-")
                    print(f"   {nm(a_)} -- {nm(b_)}: {link} | {dmin}")
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
grown_vol = 0.0; grown = np.zeros(lab.shape, bool)
seeds = np.isin(lab, acc.cand_id.astype(int).tolist()) if len(acc) else np.zeros(lab.shape, bool)
dark_p = f"{base}/work/{subj}_dark.nii.gz"
if seeds.any() and os.path.exists(dark_p):
    dark = nib.load(dark_p).get_fdata() > 0
    near = ndi.distance_transform_edt(~seeds, sampling=vox) <= GROW_MM
    grown = ndi.binary_propagation(seeds, structure=ndi.generate_binary_structure(3, 1),
                                   mask=(dark & near) | seeds)
    grown_vol = float(grown.sum() * vmm3)
    nib.save(nib.Nifti1Image(grown.astype(np.uint8), img.affine), f"{base}/review/{subj}_css_mask.nii.gz")

# v4.13 surface-based burden (van Harten et al. 2023, Discussion): blooming widens cSS PERPENDICULAR to the
# cortex far more than along it, so the cortical surface area covered - and that area as a % of the
# hemisphere's cortical surface (internally calibrated) - should depend less on field strength / TE /
# voxel size than volume does. Surface = cortex voxels facing CSF; covered = within 1.5 mm of the mask.
surf_area = {"L": 0.0, "R": 0.0}; surf_pct = {"L": 0.0, "R": 0.0}
seg_p = f"{base}/work/{subj}_seg_swispace.nii.gz"
if os.path.exists(seg_p):
    sg = nib.load(seg_p).get_fdata().astype(np.int32)
    ctx = (sg >= 1000) | np.isin(sg, [3, 42])
    tis = (sg > 0) & (sg != 24)
    surf = ctx & ndi.binary_dilation(~tis, structure=ndi.generate_binary_structure(3, 1))
    a_vox = float(np.prod(vox)) ** (2.0 / 3.0)               # mean face area of one voxel (mm2)
    cov = (ndi.distance_transform_edt(~grown, sampling=vox) <= 1.5) if grown_vol > 0 else np.zeros(sg.shape, bool)
    for h, m in (("L", ((sg >= 1000) & (sg < 2000)) | (sg == 3)), ("R", (sg >= 2000) | (sg == 42))):
        sh_ = surf & m
        surf_area[h] = float((sh_ & cov).sum() * a_vox)
        surf_pct[h] = 100.0 * (sh_ & cov).sum() / max(int(sh_.sum()), 1)

result["multifocality_0_4"] = total
result["n_foci_total"] = total_foci
n_units = total_sulci if A2 is not None else total_foci     # STRIVE-2 counts sulci
result["category"] = "absent" if n_units == 0 else ("focal" if n_units <= 3 else "disseminated")
if A2 is not None:
    result["n_sulci_total"] = total_sulci
    result["sulci"] = sulci_names
result["infratentorial_accepted"] = int(len(acc_infra))

def gyri_of(mask, hemi):
    """Destrieux gyral labels covered by a focus (nearest gyral cortex within SULC_MM, >=15 %)"""
    sl = crop(mask, SULC_MM + 6); a = A2[sl]; lo, hi = hemi_rng(hemi)
    gm = (a >= lo) & (a < hi) & ~np.isin(a, SULC) & (a % 100 != 0)
    if not gm.any(): return []
    d, ind = ndi.distance_transform_edt(~gm, sampling=vox, return_indices=True)
    m = mask[sl]; near = m & (d <= SULC_MM)
    if not near.any(): near = m
    codes = a[tuple(i[near] for i in ind)]
    cnt = np.bincount(codes - lo, minlength=100)
    return [int(c + lo) for c in np.nonzero(cnt >= max(3, 0.15 * len(codes)))[0]]

if A2 is not None:
    # Boston v2.0: two gyri belong to one focus if they touch or border the same sulcus
    gyr_adj = lambda g, t: t in touching(g) or any(is_sulc(x) for x in touching(g) & touching(t))
    n_comp = big = 0
    for h in ("L", "R"):
        G = sorted(set(g for c in acc[acc.hemi == h].cand_id.astype(int) for g in gyri_of(lab == c, h))) if len(acc) else []
        if not G: continue
        k = components(G, gyr_adj); n_comp += k
        if k == 1 and len(G) >= 3: big = 1
    result["boston2_css_lesions"] = 0 if n_comp == 0 else (1 if n_comp == 1 and not big else 2)
if A2 is not None and ICH is not None:
    result["ich_mask"] = ich_src
    result["ich_sulcal_too_close"] = ",".join(map(str, sorted(ich_excl)))
    result["ich_sulcal_check"] = ",".join(map(str, sorted(ich_check)))
# sequence (SWI vs T2*-GRE scores are not interchangeable: SWI rates higher)
# acquisition: blooming (apparent cSS size) depends on TE, field strength and the SWI processing
# (Barbosa et al., Radiol Bras 2015; van Harten et al. 2023) - record them with every score
seq = "unknown"; acq = {}
import glob as _glob
for jp in sorted(_glob.glob(f"{base}/raw/{subj}/SWI/*.json")):
    try:
        js = json.load(open(jp)); it = " ".join(js.get("ImageType", [])).upper()
        if "PHASE" in it or " P " in f" {it} ": continue
        seq = "SWI" if "SWI" in it else ("T2*-GRE" if it else seq)
        te = js.get("EchoTime"); fs = js.get("MagneticFieldStrength")
        tr = js.get("RepetitionTime"); st = js.get("SliceThickness")
        acq = dict(TE_ms=round(te * 1000, 1) if isinstance(te, (int, float)) else "",
                   TR_ms=round(tr * 1000, 1) if isinstance(tr, (int, float)) else "",
                   field_T=fs if fs is not None else "", manufacturer=js.get("Manufacturer", ""),
                   slice_thickness_mm=st if st is not None else "")
        break
    except Exception:
        pass
result["sequence"] = seq
result.update(acq)
result["voxel_mm"] = "x".join(f"{v:.2f}" for v in vox)
result["phase_available"] = int(os.path.exists(f"{base}/data/{subj}_phase.nii.gz"))
result["candidate_volume_mm3"] = round(float(acc.volume_mm3.sum()) if len(acc) else 0.0, 1)
result["grown_volume_mm3"] = round(grown_vol, 1)
for h in ("L", "R"):
    result[f"{h}_surface_mm2"] = round(surf_area[h], 1); result[f"{h}_surface_pct"] = round(surf_pct[h], 2)
result["surface_mm2"] = round(surf_area["L"] + surf_area["R"], 1)
# van Harten 2023: below ~1 mm resolution partial volume starts to bias the volume; and any size measure
# is only comparable between scans acquired with the same protocol (field, TE, voxel size)
result["volume_note"] = ("thick slices (>1.5 mm): partial volume - compare sizes only within one protocol"
                         if max(vox) > 1.5 else "compare sizes only within one protocol")
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

print(f"\n{subj}  cSS multifocality score: {total}/4   ({result['category']})   [{method}; {seq}"
      + (f" {acq['field_T']}T TE {acq['TE_ms']} ms" if acq.get("TE_ms") != "" and acq else "") + "]")
for h, nm in (("L", "left "), ("R", "right")):
    print(f"  {nm}: score {result[f'{h}_score']}  foci {result[f'{h}_foci']}  clusters {result[f'{h}_clusters']}"
          + (f"  sulci {result[f'{h}_sulci']} ({result[f'{h}_sulci_pct']}% of sulci)" if A2 is not None else ""))
if A2 is not None and sulci_names:
    print(f"  sulci: {', '.join(sulci_names)}")
print(f"  volume: candidates {result['candidate_volume_mm3']} mm3, grown (full extent) {result['grown_volume_mm3']} mm3")
if os.path.exists(seg_p):
    print(f"  cortical surface covered: {result['surface_mm2']} mm2  (L {result['L_surface_pct']} %, "
          f"R {result['R_surface_pct']} % of the hemisphere surface)  [{result['volume_note']}]")
print(f"  regions: {', '.join(result['regions'])}")
if result["n_unreviewed"]:
    print(f"  WARNING: {result['n_unreviewed']} of {result['n_candidates']} candidates were not reviewed "
          f"- score is a LOWER BOUND")
if result["infratentorial_accepted"]:
    print(f"  infratentorial siderosis (NOT in cSS score - consider classical superficial siderosis): "
          f"{result['infratentorial_accepted']} candidates")
if "boston2_css_lesions" in result:
    b2 = result["boston2_css_lesions"]
    print(f"  Boston v2.0 cSS count: " + ("0" if b2 == 0 else "1 focus (1 haemorrhagic lesion)" if b2 == 1
          else ">=2 (multifocal cSS = >=2 strictly lobar haemorrhagic lesions)"))
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

## scripts/stress_test.py

```python
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
ap.add_argument("--surface-veins", default="2", help="veins running along the cortical surface, some branching (v4.16)")
ap.add_argument("--legacy", action="store_true", help="v2 synthetic lesions (flood-filled, thicker)")
ap.add_argument("--physics", action="store_true", help="v4.18: susceptibility (dipole) model instead of darkening")
ap.add_argument("--te-ms", default="20"); ap.add_argument("--b0", default="3.0")
ap.add_argument("--calcium", default="0", help="calcification-like decoys per scan (with --physics)")
a = ap.parse_args()
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project")); S = f"{base}/scripts"
py = sys.executable
def run(*args):
    r = subprocess.run([py, *args], capture_output=True, text=True)
    if r.returncode: print(r.stdout[-500:], r.stderr[-1500:]); sys.exit("step failed: " + " ".join(args))
    return r.stdout
rows = []; feats = []; decoys = {"shown": 0, "excluded": 0, "not detected": 0}; why = {}; lost = {}
by_kind = {}                                   # decoy kind -> {"shown": n, "excluded": n, "not detected": n}
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
            run(f"{S}/make_synthetic.py", h, "8", d, s, *(["--legacy"] if a.legacy else ["--veins", a.veins, "--radial", a.radial, "--surface-veins", a.surface_veins]),
                *(["--physics", "--te-ms", a.te_ms, "--b0", a.b0, "--calcium", a.calcium] if a.physics else []))
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
                kp = f"{base}/work/{h}S_vein_kinds.csv"
                kinds = pd.read_csv(kp).set_index("vein_id")["kind"].to_dict() if os.path.exists(kp) else {}
                if os.path.exists(vp):          # which decoy kind each candidate covers most
                    vm_ = nib.load(vp).get_fdata().astype(int)
                    cdf["decoy_kind"] = [kinds.get(int(np.bincount(vm_[cm == c][vm_[cm == c] > 0]).argmax()), "")
                                         if (vm_[cm == c] > 0).any() else "" for c in cdf.cand_id]
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
                        k_ = by_kind.setdefault(kinds.get(v_, "?"), {"shown": 0, "excluded": 0, "not detected": 0})
                        if (cm[m_] > 0).any(): decoys["shown"] += 1; k_["shown"] += 1
                        elif (xm[m_] > 0).any():
                            decoys["excluded"] += 1; k_["excluded"] += 1
                            r_ = str(xr.get(int(np.bincount(xm[m_][xm[m_] > 0]).argmax()), "?")).split(":")[0].split(" (")[0]
                            why[r_] = why.get(r_, 0) + 1
                        else: decoys["not detected"] += 1; k_["not detected"] += 1
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
    for k, v in sorted(by_kind.items()):
        print(f"  {k:8s} decoys: {v['shown']} shown, {v['excluded']} excluded, {v['not detected']} not detected")
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
              "vessel_ext_mm", "vessel_wm_mm", "vessel_run_mm", "axis_normal", "local_cont_max_mm",
              "local_cont_total_mm", "terminal_ends", "local_branch_points", "score_v3", "score_v4", "score"]:
        if f in F.columns:
            L, N = F[F.is_lesion == 1][f], F[F.is_lesion == 0][f]
            V = F[F.get("is_decoy_vein", 0) == 1][f] if "is_decoy_vein" in F else N.iloc[0:0]
            Vs = F[(F.get("is_decoy_vein", 0) == 1) & (F.get("decoy_kind", "") == "surface")][f] \
                if "decoy_kind" in F else N.iloc[0:0]
            print(f"  {f:18s} cSS median {L.median():7.2f} | others median {N.median():7.2f} | AUC {auc(L, N):.2f}"
                  + (f" | vs decoy veins {auc(L, V):.2f} (n={len(V)})" if len(V) else "")
                  + (f" | vs surface veins {auc(L, Vs):.2f} (n={len(Vs)})" if len(Vs) else ""))
```

## scripts/study_meta.py

```python
#!/usr/bin/env python3
"""Study-level reference and metadata -> labels/study_meta.csv (one row per study).
usage: study_meta.py SUBJECT --status pos|neg|uncertain [--ref-L 0-2 --ref-R 0-2] [--split dev|test]
                     [--motion 0-3] [--reader NAME] [--note TEXT] [--force]
  --status   patient-level cSS reference (expert / consensus), NOT the pipeline's result
  --ref-L/R  reference multifocality per hemisphere (Charidimou 0-2)
  --split    validation split. Assigned ONCE (before the pipeline result is looked at) and never
             changed afterwards: a study that has been used to design or tune rules stays 'dev'.
             Changing an existing split needs --force and is logged in the note column.
  --motion   reader's motion / artifact grade 0 none - 3 non-diagnostic
Protocol and QC fields (vendor, field, TE, voxel size, sequence, coverage, anatomy route, QC status) are
copied from review/ID_qc.json (written by detect_css.py v4.14+); run the detector first."""
import sys, os, json, datetime, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base

if len(sys.argv) < 2 or sys.argv[1].startswith("-"): sys.exit(__doc__)
subj = sys.argv[1]; base = _base()
arg = lambda f: sys.argv[sys.argv.index(f) + 1] if f in sys.argv else None
path = f"{base}/labels/study_meta.csv"
old = pd.read_csv(path, dtype=str) if os.path.exists(path) else pd.DataFrame(columns=["study_id"])
prev = old[old.study_id == subj].iloc[0].to_dict() if (old.study_id == subj).any() else {}
row = {k: v for k, v in prev.items() if isinstance(v, str)}
row["study_id"] = subj
st = arg("--status")
if st is not None:
    if st not in ("pos", "neg", "uncertain"): sys.exit("--status must be pos, neg or uncertain")
    row["css_status"] = st
for h in ("L", "R"):
    v = arg(f"--ref-{h}")
    if v is not None:
        if v not in ("0", "1", "2"): sys.exit(f"--ref-{h} must be 0, 1 or 2")
        row[f"ref_{h}"] = v
if row.get("ref_L") and row.get("ref_R"):
    row["ref_multifocality_0_4"] = str(int(row["ref_L"]) + int(row["ref_R"]))
sp = arg("--split")
if sp is not None:
    if sp not in ("dev", "test"): sys.exit("--split must be dev or test")
    if prev.get("split") and prev["split"] != sp and isinstance(prev["split"], str):
        if "--force" not in sys.argv:
            sys.exit(f"{subj} is already in split '{prev['split']}' - the split is fixed once assigned "
                     f"(use --force only to correct a data-entry error)")
        row["note"] = (str(row.get("note", "")) + f" [split {prev['split']}->{sp} forced "
                       f"{datetime.date.today().isoformat()}]").strip()
    row["split"] = sp
    row.setdefault("split_assigned", datetime.date.today().isoformat())
for f_, k_ in (("--motion", "motion_grade"), ("--reader", "reference_reader"), ("--note", "note")):
    v = arg(f_)
    if v is not None:
        if f_ == "--motion" and v not in ("0", "1", "2", "3"): sys.exit("--motion must be 0-3")
        row[k_] = v if f_ != "--note" else (str(row.get("note", "")) + " " + v).strip()
qp = f"{base}/review/{subj}_qc.json"
if os.path.exists(qp):
    q = json.load(open(qp)); a = q.get("acquisition") or {}
    row.update(vendor=a.get("manufacturer") or "", field_T=a.get("field_T") or "", TE_ms=a.get("TE_ms") or "",
               sequence=a.get("sequence") or "", image_type=a.get("image_type") or "",
               acq_dim=a.get("acq_dim") or "", voxel_mm="x".join(f"{v:g}" for v in q["voxel_mm"]),
               slice_mm=max(q["voxel_mm"]), coverage_skull_base=int(q["coverage_includes_skull_base"]),
               anatomy=q["anatomy"], recon_all=int(q["anatomy"] == "T1 recon-all"),
               flair=int(os.path.exists(f"{base}/work/{subj}_flair_swispace.nii.gz")),
               phase=int(os.path.exists(f"{base}/data/{subj}_phase.nii.gz")),
               qc_status=q["status"], edge_agreement_mm=q.get("edge_agreement_mm"))
else:
    print(f"WARNING: {qp} not found - protocol/QC fields left empty (run detect_css.py v4.14+)")
row["updated"] = datetime.datetime.now().isoformat(timespec="seconds")
os.makedirs(os.path.dirname(path), exist_ok=True)
out = pd.concat([old[old.study_id != subj], pd.DataFrame([row])], ignore_index=True)
out.to_csv(path, index=False)
print(f"{subj}: status {row.get('css_status', '?')}, split {row.get('split', 'unassigned')}, "
      f"ref {row.get('ref_multifocality_0_4', '?')}/4, QC {row.get('qc_status', '?')} -> {path}")
```

## scripts/surface_features.py

```python
#!/usr/bin/env python3
"""Pial-surface (geodesic) features for candidates AND excluded objects - recon-all route only (v4.16).
usage: surface_features.py SUBJECT [--surface-xfm tkras_to_swi.txt] [--near-mm 3] [--reach-mm 25] [--force]
Reads subjects/ID/surf/{lh,rh}.pial, the detector output (review/ID_candidates|excluded.csv/.nii.gz) and
work/ID_dark.nii.gz; writes review/ID_surface_features.csv (one row per object: list, obj_id, fingerprint,
features) and review/ID_surface_features.json (transform, registration check). RECORDED ONLY: no rule,
score or tier uses these values until expert labels show what they are worth.

Surface -> SWI transform. FreeSurfer surfaces are in tkRAS. By default the transform is built from
subjects/ID/mri/orig.mgz (tkRAS -> T1 scanner RAS) and work/ID_swi2t1.lta (prep_anat.sh). The LTA's direction
is NOT trusted: both directions are tried and the one whose pial vertices lie on the label pial boundary of
work/ID_seg_swispace.nii.gz (which mri_vol2vol produced with the same LTA) is used. If neither agrees
within 2 mm the script stops (--force to run anyway). --surface-xfm: a 4x4 text affine tkRAS -> SWI world.

Features (each object is projected onto its nearest pial vertices, <= --near-mm away):
  surf_valid             1 if >= 3 voxels project onto the surface
  surf_normal_dist_mm    median signed distance from the pial surface along the vertex normal
                         (- inside cortex / subpial, + out in the CSF); surf_normal_iqr_mm its spread
  surf_area_mm2          pial area covered (projected vertices + 1 ring)
  surf_pieces            surface-connected pieces of the projection (>= 5 vertices); 2 when it lines both
                         banks of a sulcus (they are close in space but far apart ON the surface)
  surf_band_continuity   share of the projection in its largest surface-connected piece
  geo_length_mm          geodesic length of that piece, measured inside it (two-sweep diameter)
  geo_euclid_ratio       geodesic / straight-line distance between its two ends (>1 = folds with the cortex)
  surf_dark_cont_mm      how far DARK pial surface continues beyond the projection (<= --reach-mm);
                         vertex dark = work/ID_dark.nii.gz within 1.5 mm along its normal
  surf_exit_paths        separate dark paths leaving the projection 3-8 mm out: 0 = isolated band,
                         1 = runs on one way, 2 = runs through, >= 3 = branching
  surf_sulc / surf_curv   (v4.18) median FreeSurfer ?h.sulc (sulcal depth; > 0 = in a sulcus, < 0 = on a gyral
  surf_thickness_mm       crown) and ?h.curv (mean curvature) and ?h.thickness over the projected vertices;
  surf_sulcal_frac        share of projected vertices with sulc > 0. NaN when the file is missing. Both cSS and
                          sulcal veins sit in sulci: context for the ranking model, not a separator on its own.
  (v4.19) DIMENSIONALITY ON THE SURFACE. In one slice a vein and a cSS band both look like a line; on the pial
  surface a vein stays a thin curve (as wide as the vessel) and cSS is a patch (subpial deposit over several mm
  of cortex, across several slices). Measured on the largest surface piece:
  surf_width_mm          projected area (core vertices, no ring) / geodesic length = mean width of the footprint
  surf_cross_mm          2 x the 95th percentile geodesic distance from the footprint's centre line (the
                         diameter path), inside the piece incl. its 1-vertex ring = width at the wide parts
  surf_inradius_mm       largest geodesic distance from the edge of the piece = half the widest part; a curve
                         stays at ~1 ring (mesh spacing ~0.7-1 mm) however long it is, a patch grows
  surf_aspect            surf_cross_mm / geo_length_mm (0 = curve, ~1 = round patch)
  All NaN when the piece is shorter than 2 mm. A vein lying on a cSS patch, or cSS confined to one line of
  cortex, breaks the rule of thumb - recorded only, read with the images.
This replaces the uploaded geodesic_vein_features.py prototype: there, "continuation" was measured on the
bare mesh (always the search limit, 25 mm) and "branches" were mesh vertices with >= 6 neighbours."""
import sys, os, json, argparse, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, FP
from css_surface import apply, surface_to_swi

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("subject"); ap.add_argument("--surface-xfm", default="")
ap.add_argument("--subjects-dir", default="")
ap.add_argument("--near-mm", type=float, default=3.0); ap.add_argument("--reach-mm", type=float, default=25.0)
ap.add_argument("--force", action="store_true")
a = ap.parse_args(); subj = a.subject; base = _base()
sdir = a.subjects_dir or os.environ.get("CSS_SUBJECTS_DIR", f"{base}/subjects")
swi = nib.load(f"{base}/data/{subj}_swi.nii"); aff = swi.affine; inv_aff = np.linalg.inv(aff)
T, tsrc, checks, surfs = surface_to_swi(base, subj, sdir, a.surface_xfm, a.force)
from css_common import update_manifest, sha256
update_manifest(base, subj, "surface_features", dict(transform_source=tsrc, label_agreement_mm=checks,
    pial={h: sha256(f"{sdir}/{subj}/surf/{h}.pial") for h in ("lh", "rh")}, near_mm=a.near_mm, reach_mm=a.reach_mm))

dark = nib.load(f"{base}/work/{subj}_dark.nii.gz").get_fdata() > 0

class Surf:
    def __init__(self, V, F):
        self.V = apply(T, V); self.F = F.astype(np.int64)
        e = np.r_[F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]; e = np.unique(np.sort(e, axis=1), axis=0)
        w = np.linalg.norm(self.V[e[:, 0]] - self.V[e[:, 1]], axis=1); n = len(self.V)
        self.G = coo_matrix((np.r_[w, w], (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(n, n)).tocsr()
        fn = np.cross(self.V[F[:, 1]] - self.V[F[:, 0]], self.V[F[:, 2]] - self.V[F[:, 0]])   # area-weighted
        N = np.zeros_like(self.V)
        for k in range(3): np.add.at(N, F[:, k], fn)
        self.N = N / (np.linalg.norm(N, axis=1, keepdims=True) + 1e-9)
        if np.median(np.einsum("ij,ij->i", self.V - self.V.mean(0), self.N)) < 0: self.N = -self.N   # outward
        self.area = np.zeros(n); np.add.at(self.area, F.ravel(), np.repeat(np.linalg.norm(fn, axis=1) / 6.0, 3))
        self.tree = cKDTree(self.V)
        # vertex dark: work/ID_dark.nii.gz within 1.5 mm along the normal (subpial or just outside)
        dk = np.zeros(n, bool)
        for t in np.arange(-1.5, 1.51, 0.5):
            ijk = np.round(apply(inv_aff, self.V + t * self.N)).astype(int)
            ok = np.all((ijk >= 0) & (ijk < np.array(dark.shape)), axis=1)
            dk[ok] |= dark[tuple(ijk[ok].T)]
        self.dark = dk
S = {"L": Surf(*surfs["lh"]), "R": Surf(*surfs["rh"])}
# v4.18 per-vertex morphometry from recon-all (read only if present)
for h_, k_ in (("lh", "L"), ("rh", "R")):
    S[k_].morph = {}
    for m_ in ("sulc", "curv", "thickness"):
        p_ = f"{sdir}/{subj}/surf/{h_}.{m_}"
        if os.path.exists(p_):
            v_ = nib.freesurfer.read_morph_data(p_)
            if len(v_) == len(S[k_].V): S[k_].morph[m_] = v_

def sub_dist(sf, nodes, src, limit=np.inf):
    """Dijkstra inside the vertex subset `nodes` (bool) from `src` (indices into the full mesh)"""
    idx = np.nonzero(nodes)[0]; pos = -np.ones(len(nodes), int); pos[idx] = np.arange(len(idx))
    d = dijkstra(sf.G[idx][:, idx], indices=pos[src], limit=limit, min_only=True)
    out = np.full(len(nodes), np.inf); out[idx] = d
    return out

def sub_path(sf, nodes, a_, b_):
    """vertex indices of the shortest path a_ -> b_ inside the vertex subset `nodes`"""
    idx = np.nonzero(nodes)[0]; pos = -np.ones(len(nodes), int); pos[idx] = np.arange(len(idx))
    _, pred = dijkstra(sf.G[idx][:, idx], indices=pos[a_], return_predecessors=True)
    path, k = [], pos[b_]
    while k >= 0:
        path.append(idx[k]); k = pred[k]
    return np.array(path[::-1], int)

def features(m, hemi):
    pts = np.argwhere(m); Pw = apply(aff, pts.astype(float))
    sf = S.get(hemi) or min(S.values(), key=lambda s_: np.median(s_.tree.query(Pw)[0]))
    dd, vi = sf.tree.query(Pw)
    ok = dd <= a.near_mm
    r = dict(surf_valid=int(ok.sum() >= 3))
    if not r["surf_valid"]:
        return r
    off = np.einsum("ij,ij->i", Pw[ok] - sf.V[vi[ok]], sf.N[vi[ok]])
    r["surf_normal_dist_mm"] = round(float(np.median(off)), 2)
    r["surf_normal_iqr_mm"] = round(float(np.subtract(*np.percentile(off, [75, 25]))), 2)
    vv = np.unique(vi[ok]); mo = getattr(sf, "morph", {})
    for m_, col in (("sulc", "surf_sulc"), ("curv", "surf_curv"), ("thickness", "surf_thickness_mm")):
        r[col] = round(float(np.median(mo[m_][vv])), 3) if m_ in mo else np.nan
    r["surf_sulcal_frac"] = round(float((mo["sulc"][vv] > 0).mean()), 2) if "sulc" in mo else np.nan
    P = np.zeros(len(sf.V), bool); P[vi[ok]] = True
    P1 = P.copy(); P1[sf.G[np.nonzero(P)[0]].indices] = True          # + 1 ring: voxel grid vs mesh spacing
    r["surf_area_mm2"] = round(float(sf.area[P1].sum()), 1)
    idx = np.nonzero(P1)[0]
    nc, lab = connected_components(sf.G[idx][:, idx], directed=False)
    sizes = np.bincount(lab)
    r["surf_pieces"] = int((sizes >= 5).sum())
    r["surf_band_continuity"] = round(float(sizes.max() / len(idx)), 2)
    big = np.zeros(len(sf.V), bool); big[idx[lab == sizes.argmax()]] = True
    d0 = sub_dist(sf, big, [int(np.nonzero(big)[0][0])])
    e1 = int(np.argmax(np.where(np.isfinite(d0), d0, -1))); d1 = sub_dist(sf, big, [e1])
    e2 = int(np.argmax(np.where(np.isfinite(d1), d1, -1)))
    r["geo_length_mm"] = round(float(d1[e2]), 1)
    eu = float(np.linalg.norm(sf.V[e1] - sf.V[e2]))
    r["geo_euclid_ratio"] = round(float(d1[e2] / eu), 2) if eu > 1.0 else np.nan
    # v4.19 dimensionality on the surface: curve (vein) vs patch (cSS)
    L = float(d1[e2])
    if L >= 2.0:
        r["surf_width_mm"] = round(float(sf.area[P & big].sum() / L), 2)
        dpath = sub_dist(sf, big, sub_path(sf, big, e1, e2))
        r["surf_cross_mm"] = round(float(2 * np.percentile(dpath[big & np.isfinite(dpath)], 95)), 2)
        bi = np.nonzero(big)[0]
        edge = bi[np.array([(~big[sf.G[v].indices]).any() for v in bi])]
        if len(edge) and len(edge) < len(bi):
            de = sub_dist(sf, big, edge)
            r["surf_inradius_mm"] = round(float(de[big & np.isfinite(de)].max()), 2)
        else:
            r["surf_inradius_mm"] = np.nan
        r["surf_aspect"] = round(r["surf_cross_mm"] / L, 2)
    else:
        r.update(surf_width_mm=np.nan, surf_cross_mm=np.nan, surf_inradius_mm=np.nan, surf_aspect=np.nan)
    # dark pial surface beyond the projection, and how many separate dark paths leave it
    D = sf.dark | P1
    dc = sub_dist(sf, D, np.nonzero(P1)[0], limit=a.reach_mm)
    out = np.isfinite(dc) & ~P1
    r["surf_dark_cont_mm"] = round(float(dc[out].max()), 1) if out.any() else 0.0
    ring = out & (dc >= 3.0) & (dc <= 8.0)
    if ring.any():
        ri = np.nonzero(ring)[0]; _, rl = connected_components(sf.G[ri][:, ri], directed=False)
        r["surf_exit_paths"] = int((np.bincount(rl) >= 3).sum())
    else:
        r["surf_exit_paths"] = 0
    return r

rows = []
for lst, idc in (("candidate", "cand_id"), ("excluded", "excl_id")):
    cp, mp = f"{base}/review/{subj}_{lst if lst == 'excluded' else 'candidates'}.csv", \
        f"{base}/review/{subj}_{lst if lst == 'excluded' else 'candidates'}.nii.gz"
    if not os.path.exists(cp): continue
    d = pd.read_csv(cp)
    if not len(d): continue
    M = nib.load(mp).get_fdata().astype(int)
    for _, rr in d.iterrows():
        m = M == int(rr[idc])
        if not m.any(): continue
        rows.append(dict(list=lst, obj_id=int(rr[idc]), **{f: rr[f] for f in FP},
                         **features(m, str(rr.get("hemi", "")))))
cols = ["list", "obj_id"] + FP + ["surf_valid", "surf_normal_dist_mm", "surf_normal_iqr_mm", "surf_area_mm2",
        "surf_pieces", "surf_band_continuity", "geo_length_mm", "geo_euclid_ratio", "surf_dark_cont_mm",
        "surf_exit_paths", "surf_sulc", "surf_curv", "surf_thickness_mm", "surf_sulcal_frac",
        "surf_width_mm", "surf_cross_mm", "surf_inradius_mm", "surf_aspect"]
pd.DataFrame(rows, columns=cols).to_csv(f"{base}/review/{subj}_surface_features.csv", index=False)
json.dump(dict(subject=subj, transform_source=tsrc, label_agreement_mm=checks, tkras_to_swi=T.tolist(),
               near_mm=a.near_mm, reach_mm=a.reach_mm, note="recorded only; not validated"),
          open(f"{base}/review/{subj}_surface_features.json", "w"), indent=2)
print(f"{subj}: surface features for {len(rows)} objects -> {base}/review/{subj}_surface_features.csv")
```

## scripts/robustness.py

```python
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
```

## scripts/surface_overlay.py

```python
#!/usr/bin/env python3
"""Show candidates and excluded objects on the FreeSurfer surface, in freeview (recon-all route; v4.16).
usage: surface_overlay.py SUBJECT [--color-by list|expert] [--reader NAME] [--session N]
                          [--surface-xfm tkras_to_swi.txt] [--near-mm 3] [--force]
  --color-by list    (default) red = shown candidate, blue = excluded by a rule
  --color-by expert  the expert's calls from labels/candidate_labels.csv (expert_import.py): red cSS,
                     gold uncertain, blue surface vein (VS), cyan other veins, orange artifact,
                     grey other / normal, light green = object without a label
Every voxel of each object is projected onto its nearest pial vertex (<= --near-mm); the surface transform is
checked against the SWI labels exactly as in surface_features.py (stops above 2 mm unless --force).
Where objects overlap on the surface the SAFER class is drawn: cSS before uncertain before veins before
artifact/other - a vein can never hide a cSS call. Hidden vertices are counted per object.
Writes, in review/:
  ID_surf_<mode>.lh.annot / .rh.annot   categorical colours (annotation: no blending between categories)
  ID_surf_ids.lh.mgh / .rh.mgh          object id per vertex (+cand_id, -excl_id): click a vertex in freeview
  ID_surf_<mode>_objects.csv            per object: class, vertices, hidden vertices, Destrieux sulcus/gyrus
                                        under it (label/?h.aparc.a2009s.annot) - a check on sulcal scoring
and prints the freeview command. Display / QC only: nothing here changes candidates, calls or scores.
No model probabilities or review tiers are shown - none are validated."""
import sys, os, argparse, numpy as np, nibabel as nib, pandas as pd
from scipy.spatial import cKDTree
from nibabel.freesurfer.io import write_annot, read_annot
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, FP, match_calls, label_class
from css_surface import apply, surface_to_swi

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("subject"); ap.add_argument("--color-by", default="list", choices=["list", "expert"])
ap.add_argument("--reader", default=""); ap.add_argument("--session", type=int, default=0)
ap.add_argument("--surface-xfm", default=""); ap.add_argument("--subjects-dir", default="")
ap.add_argument("--near-mm", type=float, default=3.0); ap.add_argument("--force", action="store_true")
a = ap.parse_args(); subj = a.subject; base = _base()
sdir = a.subjects_dir or os.environ.get("CSS_SUBJECTS_DIR", f"{base}/subjects")

# classes in drawing priority (first wins where objects overlap); RGB
if a.color_by == "list":
    CLASSES = [("shown_candidate", (220, 30, 30)), ("excluded", (30, 90, 230))]
else:
    CLASSES = [("cSS", (220, 30, 30)), ("uncertain", (230, 190, 0)), ("surface_vein", (30, 90, 230)),
               ("other_vein", (0, 200, 220)), ("artifact", (240, 130, 0)), ("other_normal", (150, 150, 150)),
               ("not_labelled", (150, 230, 150))]
PRIO = {c: k for k, (c, _) in enumerate(CLASSES)}

expert_class = label_class                    # css_common (shared with consensus.py)

# ---- objects (candidates + excluded) and their class
objs = []
for lst, stem, idc in (("candidate", "candidates", "cand_id"), ("excluded", "excluded", "excl_id")):
    cp = f"{base}/review/{subj}_{stem}.csv"
    if not os.path.exists(cp): continue
    d = pd.read_csv(cp)
    if not len(d): continue
    M = nib.load(f"{base}/review/{subj}_{stem}.nii.gz").get_fdata().astype(int)
    labels = {}
    if a.color_by == "expert":
        lp = f"{base}/labels/candidate_labels.csv"
        if not os.path.exists(lp): sys.exit(f"no {lp} - import expert answers first (expert_import.py)")
        L = pd.read_csv(lp); L = L[L.study_id.astype(str) == subj]
        readers = sorted(L.reader_id.astype(str).unique())
        if not readers: sys.exit(f"{subj}: no expert labels in {lp}")
        reader = a.reader or (readers[0] if len(readers) == 1 else sys.exit(f"several readers {readers}: --reader"))
        if reader not in readers: sys.exit(f"{subj}: no labels from reader '{reader}' (have {readers})")
        L = L[L.reader_id.astype(str) == reader]
        ses = a.session or int(L.session.max()); L = L[(L.session == ses) & (L.list == lst)]
        if len(L):
            fp = L.cand_uid.str.split(":", expand=True).iloc[:, 1:5].astype(float)
            fp.columns = FP
            m, _ = match_calls(d.rename(columns={idc: "cand_id"}), fp.assign(call=L.label.values))
            labels = m or {}
    for _, r in d.iterrows():
        i = int(r[idc])
        cls = ("shown_candidate" if lst == "candidate" else "excluded") if a.color_by == "list" else \
            (expert_class(labels[i]) if i in labels else "not_labelled")
        objs.append(dict(list=lst, obj_id=i, sid=i if lst == "candidate" else -i, hemi=str(r.get("hemi", "")),
                         cls=cls, label=labels.get(i, ""), mask=(M == i), **{f: r[f] for f in FP}))
if not objs: sys.exit(f"{subj}: no candidates or excluded objects")

# ---- surfaces in SWI space (checked transform), projection, safer-class-wins drawing
T, tsrc, checks, surfs = surface_to_swi(base, subj, sdir, a.surface_xfm, a.force)
swi = nib.load(f"{base}/data/{subj}_swi.nii"); aff = swi.affine
out = {}
for h, H in (("lh", "L"), ("rh", "R")):
    V, F = surfs[h]; Vw = apply(T, V); tree = cKDTree(Vw)
    nb = [[] for _ in range(len(V))]
    for e0, e1 in np.r_[F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]:
        nb[e0].append(e1); nb[e1].append(e0)
    a2p = f"{sdir}/{subj}/label/{h}.aparc.a2009s.annot"
    a2 = read_annot(a2p) if os.path.exists(a2p) else None
    out[h] = dict(n=len(V), best=np.full(len(V), len(CLASSES)), sid=np.zeros(len(V), np.float32), tree=tree,
                  nb=nb, a2=a2)
for o in objs:
    pts = apply(aff, np.argwhere(o["mask"]).astype(float))
    hs = ["lh"] if o["hemi"] == "L" else ["rh"] if o["hemi"] == "R" else ["lh", "rh"]
    best_h = min(hs, key=lambda h_: np.median(out[h_]["tree"].query(pts)[0]))
    dd, vi = out[best_h]["tree"].query(pts)
    core = np.unique(vi[dd <= a.near_mm])
    ring = np.unique(np.r_[core, [n for v in core for n in out[best_h]["nb"][v]]]).astype(int) if len(core) else core
    o.update(hemi_surf=best_h, core=core, ring=ring, mean_dist_mm=round(float(np.median(dd)), 2))
    # draw: a vertex keeps the safest (lowest-priority-number) class claiming it
    S = out[best_h]; p = PRIO[o["cls"]]; win = p < S["best"][ring]
    S["best"][ring[win]] = p; S["sid"][ring[win]] = o["sid"]
rows = []
for o in objs:
    S = out[o["hemi_surf"]]
    sulc = gyr = ""
    if S["a2"] is not None and len(o["core"]):
        names = np.array([n.decode() if isinstance(n, bytes) else n for n in S["a2"][2]])
        li = S["a2"][0][o["core"]]; nm = names[li[li >= 0]].astype(str)        # -1 = vertex without a label
        s_ = nm[np.char.startswith(nm, "S_") | (np.char.find(nm, "Lat_Fis") >= 0)]
        g_ = nm[np.char.startswith(nm.astype(str), "G_")]
        sulc = pd.Series(s_).mode().iloc[0] if len(s_) else ""; gyr = pd.Series(g_).mode().iloc[0] if len(g_) else ""
    rows.append(dict(list=o["list"], obj_id=o["obj_id"], **{f: o[f] for f in FP}, hemi=o["hemi_surf"][0].upper(),
                     cls=o["cls"], expert_label=o["label"], n_vertices=len(o["core"]),
                     n_hidden=int((S["best"][o["core"]] < PRIO[o["cls"]]).sum()) if len(o["core"]) else 0,
                     median_dist_to_pial_mm=o["mean_dist_mm"], sulcus=sulc, gyrus=gyr))
# ---- write annotation (categorical), id overlay, table; print the freeview command
ctab = np.array([[*rgb, 0, 0] for _, rgb in CLASSES], np.int32)
ctab[:, 4] = ctab[:, 0] + ctab[:, 1] * 256 + ctab[:, 2] * 65536
mode = a.color_by; files = {}
for h, S in out.items():
    lab = np.where(S["best"] < len(CLASSES), S["best"], -1).astype(np.int32)
    ap_ = f"{base}/review/{subj}_surf_{mode}.{h}.annot"; ip_ = f"{base}/review/{subj}_surf_ids.{h}.mgh"
    write_annot(ap_, lab, ctab, [c for c, _ in CLASSES], fill_ctab=False)
    nib.save(nib.MGHImage(S["sid"].reshape(-1, 1, 1), np.eye(4)), ip_)
    files[h] = (ap_, ip_)
tab = pd.DataFrame(rows); tab.to_csv(f"{base}/review/{subj}_surf_{mode}_objects.csv", index=False)
hid = tab[(tab.n_hidden > 0)]
print(f"{subj}: {len(tab)} objects on the surface ({', '.join(f'{v} {k}' for k, v in tab.cls.value_counts().items())})"
      f"; {int((tab.n_vertices == 0).sum())} too far from the pial surface (> {a.near_mm:g} mm)"
      + (f"; {len(hid)} partly hidden under a safer class" if len(hid) else ""))
sv = "inflated" if os.path.exists(f"{sdir}/{subj}/surf/lh.inflated") else "pial"
print("view:  freeview -viewport 3d " + " ".join(
    f"-f {sdir}/{subj}/surf/{h}.{sv}:annot={files[h][0]}:overlay={files[h][1]}:overlay_threshold=0.5,999"
    for h in ("lh", "rh")))
print("       (annotation = colours; the overlay holds object ids: + candidate, - excluded; click a vertex to read it)")
```

## tests/check_v418.py

```python
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
# 7. robustness.py: an identity perturbation (gamma 1) must reproduce the copy exactly; copies are cleaned up
import subprocess, glob
r = subprocess.run([sys.executable, f"{B}/scripts/robustness.py", "PH3", "--perturb", "copy,gamma:1.0,shift:0",
                    "--old", f"{B}/work/PH3_truth.nii.gz", "--ids", "1,2,3,4,5,6"], capture_output=True, text=True)
if r.returncode: print(r.stdout[-800:], r.stderr[-1500:]); fail("robustness.py failed")
rb = pd.read_csv(f"{B}/review/PH3_robustness.csv").set_index("perturbation")
g, c0 = rb.loc["gamma:1"], rb.loc["copy"]
if not (g.candidates == c0.candidates == g.found_again and g.new == 0 and g.spearman == 1.0):
    fail(f"identity perturbation changed the result: {g.to_dict()}")
if glob.glob(f"{B}/*/PH3_R*"): fail("robustness copies not cleaned up")
print(f"   robustness: copy {int(c0.candidates)} candidates, gamma 1 identical; 1-voxel shift -> {int(rb.loc['shift:0'].candidates)} "
      f"candidates, {int(rb.loc['shift:0'].found_again)}/{int(c0.candidates)} found again (phantom)")
print("RESULT v4.18 OK")
```

## scripts/make_swi.py

```python
#!/usr/bin/env python3
"""Build an SWI image from gradient-echo MAGNITUDE + PHASE (v4.18) - for datasets that store the raw GRE
(e.g. OpenNeuro ds000221 / MPI-Leipzig LEMON: acq-mag_GRE + acq-phase_GRE) instead of a processed SWI.
usage: make_swi.py MAG.nii.gz PHASE.nii.gz OUT_swi.nii.gz [--sigma-mm 4] [--power 4] [--phase-sign auto|+1|-1]
                   [--hp-out OUT_hpphase.nii.gz]
Standard SWI processing (Haacke 2004): phase -> radians, homodyne high-pass per axial slice
(angle(z * conj(lowpass z)), Gaussian low-pass sigma --sigma-mm, wrap-safe; the same filter the reviewer
uses), a phase mask that darkens PARAMAGNETIC voxels (veins, haemosiderin) linearly to 0 at pi, and
SWI = magnitude x mask^power. Nothing else is changed: the output keeps the magnitude's grid and affine.

Phase sign: which sign is paramagnetic differs between vendors and reconstructions. With --phase-sign auto
(default) it is read from the image: paramagnetic voxels are also dark on the magnitude at this TE, so the
sign of the median high-pass phase of the darkest-for-their-neighbourhood brain voxels is used. The
evidence is printed; if it is weak the script stops and asks for --phase-sign. Check the result once per
protocol: veins must be DARKER on the SWI than on the magnitude (calcium is not darkened).
Research software, not for clinical reading."""
import sys, os, argparse, numpy as np, nibabel as nib
from scipy import ndimage as ndi

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("mag"); ap.add_argument("phase"); ap.add_argument("out")
ap.add_argument("--sigma-mm", type=float, default=4.0); ap.add_argument("--power", type=int, default=4)
ap.add_argument("--phase-sign", default="auto"); ap.add_argument("--hp-out", default="")
a = ap.parse_args()

mi, pi_ = nib.load(a.mag), nib.load(a.phase)
if mi.shape[:3] != pi_.shape[:3] or not np.allclose(mi.affine, pi_.affine, atol=1e-3):
    sys.exit("magnitude and phase are not on the same grid - they must come from the same acquisition")
M = np.asarray(mi.dataobj, np.float32); P = np.asarray(pi_.dataobj, np.float32)
if M.ndim == 4 or P.ndim == 4:
    sys.exit("4-D (multi-echo) input: pass one echo (e.g. the longest usable TE) as 3-D images")
vox = np.array(mi.header.get_zooms()[:3], float)
si = int(np.argmax(np.abs(mi.affine[2, :3])))          # superior-inferior array axis: slice-wise filter
lo, hi = float(np.nanmin(P)), float(np.nanmax(P))
if hi - lo > 2 * np.pi + 0.5:                            # scanner units (e.g. -4096..4095 or 0..4095)
    P = (P - lo) / (hi - lo) * 2 * np.pi - np.pi
    unit = f"scanner units {lo:g}..{hi:g} -> radians"
else:
    unit = "radians"

# homodyne high-pass, slice by slice in the plane perpendicular to the S-I axis
z = np.exp(1j * P.astype(np.float64))
inpl = [ax for ax in range(3) if ax != si]
sig = [0.0, 0.0, 0.0]
for ax in inpl: sig[ax] = a.sigma_mm / vox[ax]
lp = ndi.gaussian_filter(z.real, sig) + 1j * ndi.gaussian_filter(z.imag, sig)     # sigma 0 on the S-I axis
HP = np.angle(z * np.conj(lp)).astype(np.float32)

# brain-ish mask for the sign estimate: bright magnitude, eroded away from the head edge
thr = 0.15 * float(np.percentile(M[M > 0], 99)) if (M > 0).any() else 0.0
head = ndi.binary_fill_holes(M > thr)
lab, n = ndi.label(head)
if n: head = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
core = ndi.binary_erosion(head, iterations=max(1, int(round(10 / vox.min()))))
if a.phase_sign == "auto":
    rel = M / (ndi.gaussian_filter(M, 3.0 / vox) + 1e-6)
    dark = core & (rel < 0.75) & (np.abs(HP) > 0.2)
    if dark.sum() < 200:
        sys.exit(f"too few dark brain voxels ({int(dark.sum())}) to read the phase sign - give --phase-sign +1 or -1")
    med = float(np.median(HP[dark])); frac = float((np.sign(HP[dark]) == np.sign(med)).mean())
    if abs(med) < 0.05 or frac < 0.6:
        sys.exit(f"phase sign unclear (median {med:+.2f} rad, {frac:.0%} agree) - give --phase-sign +1 or -1")
    sign = 1 if med > 0 else -1
    how = f"auto: dark-on-magnitude voxels median high-pass phase {med:+.2f} rad, {frac:.0%} same sign, n={int(dark.sum())}"
else:
    sign = int(a.phase_sign); how = "given"
    if sign not in (1, -1): sys.exit("--phase-sign must be auto, +1 or -1")

# phase mask: paramagnetic side (sign) goes linearly from 1 at 0 to 0 at pi
para = sign * HP
mask = np.where(para > 0, np.clip(1.0 - para / np.pi, 0.0, 1.0), 1.0).astype(np.float32)
S = (M * mask ** a.power).astype(np.float32)
out = nib.Nifti1Image(S, mi.affine, mi.header); out.set_data_dtype(np.float32)
nib.save(out, a.out)
if a.hp_out:
    nib.save(nib.Nifti1Image(HP, mi.affine), a.hp_out)
red = float(1 - np.median(S[core & (mask < 0.9)] / (M[core & (mask < 0.9)] + 1e-6))) if (core & (mask < 0.9)).any() else 0.0
print(f"{os.path.basename(a.out)}: SWI = magnitude x phase-mask^{a.power}; phase {unit}; high-pass sigma {a.sigma_mm} mm "
      f"in-plane (S-I axis {si}); paramagnetic sign {sign:+d} ({how}); {(core & (mask < 0.9)).mean() if core.any() else 0:.1%} "
      f"of brain voxels darkened, median by {red:.0%}")
```

## scripts/run_lemon.sh

```bash
#!/bin/zsh
# Batch for OpenNeuro ds000221 (MPI-Leipzig LEMON): healthy volunteers with gradient-echo magnitude + phase,
# MP2RAGE T1 and 3-D FLAIR. Use: cSS-NEGATIVE CONTROL scans (false positives per scan, robustness) - 'dev' only.
# usage: run_lemon.sh LEMON_DIR [--download] [--parallel N] [--subjects "sub-010239 sub-010255"] [--summary]
#   LEMON_DIR   folder with sub-XXXXXX/ses-01/anat/ (created by --download); no spaces in the path
#   --download  fetch the 4 files per subject from the public OpenNeuro S3 bucket (curl; skips files present)
#   --parallel  subjects processed at the same time (default 1 since v4.19.1; 2 crashed a laptop)
#   --threads   CPU threads per subject for recon-all / SynthSeg (default 4; CSS_THREADS)
#   --subjects  default: the 26 volunteers aged 55-75 with GRE + 3-D FLAIR
#   --summary   only print the summary table of finished subjects
#   --force     re-run subjects that already finished (default: skipped, so an interrupted batch can be restarted)
# Per subject: make_swi.py (magnitude + phase -> SWI); mp2rage_robust.py (UNI + INV1 + INV2 -> T1 with a clean
#   background: the noisy MP2RAGE background made recon-all fail on P201) -> run_css.sh ID swi --t1 --flair --recon
#   --phase; if recon-all still fails: the same without --recon (SynthSeg on the T1 + mri_coreg; no surfaces)
#   -> study_meta.py ID --split dev. IDs P201... are fixed in labels/id_map.csv (study code <-> source id).
# Logs: logs/ID_run.log. Summary: review/lemon_summary.csv.
# These volunteers are PRESUMED cSS-negative (community volunteers, not expert-read) - never claim more.
B=${CSS_BASE:-$HOME/css_project}; SC=$B/scripts
DEFAULT="sub-010239 sub-010255 sub-010257 sub-010258 sub-010146 sub-010148 sub-010152 sub-010267 sub-010268
sub-010272 sub-010275 sub-010276 sub-010277 sub-010278 sub-010279 sub-010280 sub-010282 sub-010283 sub-010284
sub-010285 sub-010286 sub-010287 sub-010288 sub-010289 sub-010290 sub-010291"
FILES="acq-mag_GRE acq-phase_GRE acq-mp2rage_T1w inv-1_mp2rage inv-2_mp2rage acq-highres_FLAIR"
URL=https://s3.amazonaws.com/openneuro.org/ds000221

if [ "$1" = "--one" ]; then                      # internal: one subject (run in parallel by xargs)
  s=$2; ID=$3; p=$LEMON_DIR/$s/ses-01/anat/${s}_ses-01
  mkdir -p $B/work $B/logs
  if [ -z "$LEMON_FORCE" ] && grep -q "^== $ID done" $B/logs/${ID}_run.log 2>/dev/null; then
    echo "SKIP   $ID ($s) - already done (--force to re-run)"; exit 0
  fi
  if ( echo "== $ID = $s ($(date))"
       for f in $(echo $FILES); do [ -s ${p}_${f}.nii.gz ] || { echo "ERROR: missing ${p}_${f}.nii.gz"; exit 1; }; done
       python $SC/make_swi.py ${p}_acq-mag_GRE.nii.gz ${p}_acq-phase_GRE.nii.gz $B/work/${ID}_gre_swi.nii.gz || exit 1
       python $SC/mp2rage_robust.py ${p}_acq-mp2rage_T1w.nii.gz ${p}_inv-1_mp2rage.nii.gz ${p}_inv-2_mp2rage.nii.gz \
         $B/work/${ID}_t1_robust.nii.gz || exit 1
       ANAT="recon-all"
       if ! NOVIEW=1 zsh $SC/run_css.sh $ID $B/work/${ID}_gre_swi.nii.gz --t1 $B/work/${ID}_t1_robust.nii.gz \
              --flair ${p}_acq-highres_FLAIR.nii.gz --recon --phase ${p}_acq-phase_GRE.nii.gz; then
         echo "== recon-all route failed for $ID - falling back to SynthSeg on the T1 (no surfaces, mri_coreg)"
         [ -d $B/subjects/$ID ] && mv $B/subjects/$ID $B/subjects/${ID}_failed_$(date +%Y%m%d_%H%M%S)
         ANAT="SynthSeg fallback (recon-all failed)"
         NOVIEW=1 zsh $SC/run_css.sh $ID $B/work/${ID}_gre_swi.nii.gz --t1 $B/work/${ID}_t1_robust.nii.gz \
           --flair ${p}_acq-highres_FLAIR.nii.gz --phase ${p}_acq-phase_GRE.nii.gz || exit 1
       fi
       python $SC/study_meta.py $ID --split dev \
         --note "LEMON ds000221 $s healthy volunteer: PRESUMED cSS-negative, not expert-read; GRE->SWI by make_swi.py; T1 robust MP2RAGE; anatomy $ANAT" || exit 1
       echo "== $ID done ($(date))" ) > $B/logs/${ID}_run.log 2>&1; then
    echo "OK     $ID ($s)"
  else
    echo "FAILED $ID ($s) - see logs/${ID}_run.log"; tail -3 $B/logs/${ID}_run.log
  fi
  exit 0
fi

L=$1; shift 2>/dev/null || true
[ -n "$L" ] || { sed -n 2,14p $0; exit 1; }
DL=0; PAR=1; SUBS=$DEFAULT; SUMONLY=0; export CSS_THREADS=${CSS_THREADS:-4}
while [ $# -gt 0 ]; do
  case "$1" in
    --download) DL=1; shift;;
    --parallel) PAR=$2; shift 2;;
    --threads) export CSS_THREADS=$2; shift 2;;
    --subjects) SUBS=$2; shift 2;;
    --summary) SUMONLY=1; shift;;
    --force) export LEMON_FORCE=1; shift;;
    *) echo "run_lemon.sh: unknown option $1"; exit 1;;
  esac
done
case "$L" in *" "*) echo "ERROR: no spaces in LEMON_DIR please ($L)"; exit 1;; esac
mkdir -p $L $B/labels $B/logs $B/review
export LEMON_DIR=$(cd $L && pwd)

# fixed study codes: reuse labels/id_map.csv, append new subjects as the next P2xx
MAP=$B/labels/id_map.csv
[ -f $MAP ] || echo "study_id,source_dataset,source_id" > $MAP
for s in $(echo $SUBS); do
  if ! grep -q ",ds000221,$s\$" $MAP; then
    n=$(grep -c ",ds000221," $MAP); printf "P%03d,ds000221,%s\n" $((201 + n)) $s >> $MAP
  fi
done

if [ $SUMONLY -eq 0 ]; then
  if [ $DL -eq 1 ]; then
    for s in $(echo $SUBS); do for f in $(echo $FILES); do
      o=$LEMON_DIR/$s/ses-01/anat/${s}_ses-01_${f}.nii.gz; mkdir -p $(dirname $o)
      [ -s $o ] || curl -fsSL -o $o $URL/$s/ses-01/anat/${s}_ses-01_${f}.nii.gz || { echo "download FAILED: $s $f"; rm -f $o; }
    done; done
  fi
  echo "LEMON: $(echo $SUBS | wc -w | tr -d ' ') subjects, $PAR at a time, $CSS_THREADS threads, low priority (nice); recon-all ~1-3 h each (progress: tail -f $B/logs/<ID>_run.log)"
  for s in $(echo $SUBS); do echo "$s $(grep ",ds000221,$s\$" $MAP | cut -d, -f1)"; done | xargs -n 2 -P $PAR nice -n 15 zsh $0 --one
fi

python - "$B" $(echo $SUBS) <<'EOF'
import sys, os, json, pandas as pd
B = sys.argv[1]; subs = sys.argv[2:]
m = pd.read_csv(f"{B}/labels/id_map.csv"); m = m[m.source_dataset == "ds000221"].set_index("source_id")
rows = []
for s in subs:
    if s not in m.index: continue
    i = m.loc[s, "study_id"]; r = dict(study_id=i, source_id=s)
    c, x, q = (f"{B}/review/{i}_{k}" for k in ("candidates.csv", "excluded.csv", "qc.json"))
    if not os.path.exists(c): r["status"] = "not run / failed"; rows.append(r); continue
    cd = pd.read_csv(c); r.update(status="done", candidates=len(cd), excluded=len(pd.read_csv(x)) if os.path.exists(x) else None)
    if "review_tier" in cd: r.update(tier1=int((cd.review_tier == 1).sum()), tier2=int((cd.review_tier == 2).sum()))
    if os.path.exists(q):
        qq = json.load(open(q)); r.update(qc=qq.get("status"), label_shift_mm=qq.get("label_shift_mm"), anatomy=qq.get("anatomy"))
    rows.append(r)
d = pd.DataFrame(rows)
if len(d):
    d.to_csv(f"{B}/review/lemon_summary.csv", index=False)
    print(d.to_string(index=False))
    ok = d[d.status == "done"]
    if len(ok):
        print(f"\ncandidates per PRESUMED-negative scan: median {ok.candidates.median():g}, range {ok.candidates.min():g}-"
              f"{ok.candidates.max():g} (n={len(ok)}); QC ok in {int((ok.get('qc') == 'ok').sum())}/{len(ok)}. They count as "
              f"false positives only if an expert confirms the scan is cSS-free.")
EOF
```

## tests/check_make_swi.py

```python
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
```

## scripts/evaluate.py

```python
#!/usr/bin/env python3
"""Patient-grouped evaluation of the detector against expert labels (v4.18; protocol after VALDO / CenSynCMB 2026:
per-patient counts, patient bootstrap, stratified by protocol). Reads labels/, never changes them.
usage: evaluate.py [--split test|dev|all] [--source consensus|reader:NAME] [--k 5,10,20] [--boot 1000] [--seed 1]
                   [--by field_T,sequence,slice_group] [--tag NAME]
Reference standard, per study (labels/study_meta.csv gives split and patient status):
  cSS lesions = objects called cSS (labels/consensus.csv, or one reader's latest session in candidate_labels.csv)
                among shown candidates AND excluded objects (sheets made with --include-excluded),
                + reader-outlined missed cSS (labels/missed_lesions.csv; outlines within 5 mm across readers count once)
  detected    = the lesion is a SHOWN candidate in the current detector output (matched by cand_uid fingerprint);
                its rank = cand_id. Missed lesions are split by why: excluded (rule id) / dropped / not generated.
  false positive = shown candidate called anything but cSS or uncertain; 'uncertain' calls are counted apart.
Reported (with 95 % patient-bootstrap intervals): lesion sensitivity (all shown and top-k), patient sensitivity
(cSS-positive patients with >= 1 lesion shown), false positives per scan and per cSS-NEGATIVE scan, uncertain per
scan; a FROC table (rank cut-off -> sensitivity vs mean FP per scan); a per-study table; strata from --by.
Labels that no longer match the current detector output (re-run after labelling) are reported, not guessed.
Writes review/evaluation_<tag>.json, _patients.csv, _froc.csv. Default split is TEST: dev studies were used to
design the rules, so their numbers are optimistic by construction."""
import sys, os, json, argparse, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, cand_uid, label_class

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--split", default="test", choices=["test", "dev", "all"])
ap.add_argument("--source", default="consensus"); ap.add_argument("--k", default="5,10,20")
ap.add_argument("--boot", type=int, default=1000); ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--by", default=""); ap.add_argument("--tag", default="")
a = ap.parse_args(); B = _base(); L = f"{B}/labels"
KS = [int(x) for x in a.k.split(",") if x]
tag = a.tag or f"{a.split}_{a.source.replace(':', '-')}"

meta = pd.read_csv(f"{L}/study_meta.csv", dtype=str)
if a.split != "all": meta = meta[meta.split == a.split]
if not len(meta): sys.exit(f"no studies with split '{a.split}' in labels/study_meta.csv")

if a.source == "consensus":
    lab = pd.read_csv(f"{L}/consensus.csv")
    pend = lab[lab.method == "pending"]; lab = lab[lab.method != "pending"]
    lab = lab.assign(cls=lab.consensus_class)
elif a.source.startswith("reader:"):
    rd = a.source.split(":", 1)[1]; lab = pd.read_csv(f"{L}/candidate_labels.csv")
    lab = lab[lab.reader_id == rd]
    if not len(lab): sys.exit(f"no labels from reader {rd}")
    lab = lab[lab.session == lab.groupby("study_id").session.transform("max")]
    pend = lab.iloc[0:0]; lab = lab.assign(cls=lab.label.map(label_class))
else:
    sys.exit("--source consensus or reader:NAME")
mp = f"{L}/missed_lesions.csv"
missed = pd.read_csv(mp) if os.path.exists(mp) else pd.DataFrame(columns=["study_id"])

rows, warn = [], []
for _, st in meta.iterrows():
    s = st.study_id; cp, xp = f"{B}/review/{s}_candidates.csv", f"{B}/review/{s}_excluded.csv"
    if not os.path.exists(cp): warn.append(f"{s}: no detector output - skipped"); continue
    c = pd.read_csv(cp); x = pd.read_csv(xp) if os.path.exists(xp) else pd.DataFrame(columns=c.columns)
    uid_rank = {cand_uid(s, r): int(r.cand_id) for _, r in c.iterrows()}
    uid_rule = {cand_uid(s, r): str(r.get("rule_version", r.get("reason", ""))) for _, r in x.iterrows()}
    ls = lab[lab.study_id == s]
    lesions, fp_ranks, unc_ranks, stale = [], [], [], 0
    seen_shown = set()
    for _, r in ls.iterrows():
        u = r.cand_uid
        if u in uid_rank:
            rk = uid_rank[u]; seen_shown.add(rk)
            if r.cls == "cSS": lesions.append(dict(fate="shown", rank=rk))
            elif r.cls == "uncertain": unc_ranks.append(rk)
            else: fp_ranks.append(rk)
        elif u in uid_rule:
            if r.cls == "cSS": lesions.append(dict(fate="excluded", rank=None, why=uid_rule[u]))
        else:
            stale += 1
    unlabelled = sorted(set(uid_rank.values()) - seen_shown)
    # reader-outlined missed lesions (fates from label_missed.py), outlines from several readers counted once
    ms = missed[missed.study_id == s] if len(missed) else missed
    kept = []
    for _, r in ms.iterrows():
        p = np.array([r.vox_i, r.vox_j, r.vox_k], float)
        if any(np.linalg.norm(p - q) <= 5.0 for q in kept): continue
        kept.append(p)
        f = str(r.fate)
        if f == "shown":
            ob = int(r.obj_id) if pd.notna(r.obj_id) else None
            if ob is not None and any(l_["rank"] == ob for l_ in lesions): continue    # already counted
            lesions.append(dict(fate="shown", rank=ob))
            if ob in fp_ranks: fp_ranks.remove(ob)                 # it overlaps real cSS: not a false positive
        elif f == "excluded": lesions.append(dict(fate="excluded", rank=None, why=str(r.rule_version)))
        else: lesions.append(dict(fate=f, rank=None, why=str(r.reason)))
    row = dict(study_id=s, split=st.split, css_status=st.get("css_status"), n_candidates=len(c),
               n_ref_lesions=len(lesions), n_detected=sum(l_["fate"] == "shown" for l_ in lesions),
               n_lost_excluded=sum(l_["fate"] == "excluded" for l_ in lesions),
               n_lost_dropped=sum(l_["fate"] == "dropped" for l_ in lesions),
               n_lost_not_generated=sum(l_["fate"] == "not generated" for l_ in lesions),
               n_fp=len(fp_ranks), n_uncertain=len(unc_ranks), n_unlabelled=len(unlabelled), n_stale_labels=stale,
               n_pending=int((pend.study_id == s).sum()) if len(pend) else 0,
               lesion_ranks=";".join(str(l_["rank"]) for l_ in lesions if l_["fate"] == "shown"),
               fp_ranks=";".join(map(str, sorted(fp_ranks))),
               lost=";".join(f"{l_['fate']}:{l_.get('why', '')}" for l_ in lesions if l_["fate"] != "shown"))
    for k in KS:
        row[f"det_top{k}"] = sum(1 for l_ in lesions if l_["fate"] == "shown" and l_["rank"] is not None and l_["rank"] <= k)
        row[f"fp_top{k}"] = sum(1 for rk in fp_ranks if rk <= k)
    for col in ("field_T", "sequence", "vendor", "voxel_mm", "slice_mm", "qc_status", "anatomy"):
        row[col] = st.get(col)
    try: sm = float(st.get("slice_mm"))
    except (TypeError, ValueError): sm = np.nan
    row["slice_group"] = "unknown" if np.isnan(sm) else ("<=1.2 mm" if sm <= 1.2 else ("<=2 mm" if sm <= 2.0 else ">2 mm"))
    if stale: warn.append(f"{s}: {stale} label(s) match no current candidate/excluded object (detector re-run after labelling?)")
    if unlabelled: warn.append(f"{s}: {len(unlabelled)} shown candidate(s) without a label - neither TP nor FP")
    rows.append(row)
P = pd.DataFrame(rows)
if not len(P): sys.exit("no evaluable studies\n" + "\n".join(warn))

def metrics(d):
    ref = d.n_ref_lesions.sum(); pos = d[d.n_ref_lesions > 0]; neg = d[d.css_status == "neg"]
    m = dict(n_studies=len(d), n_ref_lesions=int(ref),
             lesion_sens=d.n_detected.sum() / ref if ref else np.nan,
             patient_sens=float((pos.n_detected > 0).mean()) if len(pos) else np.nan,
             fp_per_scan=float(d.n_fp.mean()), fp_per_neg_scan=float(neg.n_fp.mean()) if len(neg) else np.nan,
             uncertain_per_scan=float(d.n_uncertain.mean()))
    for k in KS:
        m[f"sens_top{k}"] = d[f"det_top{k}"].sum() / ref if ref else np.nan
        m[f"fp_per_scan_top{k}"] = float(d[f"fp_top{k}"].mean())
    return m

def boot(d, rng):
    out = {}
    if len(d) < 2: return out
    bs = [metrics(d.iloc[rng.integers(0, len(d), len(d))]) for _ in range(a.boot)]
    for key in bs[0]:
        v = np.array([b_[key] for b_ in bs], float); v = v[np.isfinite(v)]
        if len(v) and key not in ("n_studies", "n_ref_lesions"):
            out[key] = [round(float(np.percentile(v, 2.5)), 3), round(float(np.percentile(v, 97.5)), 3)]
    return out

rng = np.random.default_rng(a.seed)
overall = metrics(P); ci = boot(P, rng)
# FROC: rank cut-off k -> pooled lesion sensitivity vs mean FP per scan (ranks are per study)
maxk = int(P.n_candidates.max())
froc = []
for k in range(1, maxk + 1):
    det = sum(sum(1 for t in str(r.lesion_ranks).split(";") if t not in ("", "None") and int(t) <= k) for _, r in P.iterrows())
    fp = np.mean([sum(1 for t in str(r.fp_ranks).split(";") if t and int(t) <= k) for _, r in P.iterrows()])
    froc.append(dict(top_k=k, lesion_sens=det / overall["n_ref_lesions"] if overall["n_ref_lesions"] else np.nan,
                     fp_per_scan=float(fp)))
strata = {}
for col in [c_ for c_ in a.by.split(",") if c_]:
    if col not in P: warn.append(f"--by {col}: not a column"); continue
    strata[col] = {str(g): dict(metrics(d), ci=boot(d, rng)) for g, d in P.groupby(P[col].astype(str))}

os.makedirs(f"{B}/review", exist_ok=True)
P.to_csv(f"{B}/review/evaluation_{tag}_patients.csv", index=False)
pd.DataFrame(froc).to_csv(f"{B}/review/evaluation_{tag}_froc.csv", index=False)
json.dump(dict(split=a.split, source=a.source, overall=overall, ci95=ci, strata=strata, warnings=warn,
               pending_adjudication=int(len(pend)), boot=a.boot, seed=a.seed,
               note="patient-grouped; research software; dev-split numbers are optimistic (rules were tuned on them)"),
          open(f"{B}/review/evaluation_{tag}.json", "w"), indent=2, default=float)

f = lambda v: "n/a" if v is None or (isinstance(v, float) and not np.isfinite(v)) else (f"{v:.2f}" if isinstance(v, float) else str(v))
fci = lambda k_: f" [{ci[k_][0]:.2f}-{ci[k_][1]:.2f}]" if k_ in ci else ""
print(f"evaluation ({a.split} split, reference: {a.source}): {overall['n_studies']} studies, {overall['n_ref_lesions']} reference cSS lesions")
if a.split != "test": print("  NOTE: not the test split - these numbers are optimistic (rules were designed on dev studies)")
if overall["n_studies"] < 10: print(f"  NOTE: only {overall['n_studies']} studies - intervals are wide and unstable")
print(f"  lesion sensitivity (all shown)  {f(overall['lesion_sens'])}{fci('lesion_sens')}")
for k in KS:
    print(f"  lesion sensitivity, top {k:<3d}     {f(overall[f'sens_top{k}'])}{fci(f'sens_top{k}')}   FP per scan in top {k}: "
          f"{f(overall[f'fp_per_scan_top{k}'])}{fci(f'fp_per_scan_top{k}')}")
print(f"  patient sensitivity             {f(overall['patient_sens'])}{fci('patient_sens')}")
print(f"  false positives per scan        {f(overall['fp_per_scan'])}{fci('fp_per_scan')}   per cSS-negative scan "
      f"{f(overall['fp_per_neg_scan'])}{fci('fp_per_neg_scan')}   uncertain per scan {f(overall['uncertain_per_scan'])}")
lost = P[["n_lost_excluded", "n_lost_dropped", "n_lost_not_generated"]].sum()
print(f"  lesions not shown: {int(lost.n_lost_excluded)} excluded by a rule, {int(lost.n_lost_dropped)} dropped at a gate, "
      f"{int(lost.n_lost_not_generated)} never generated (details: lost column of the patient table)")
for col, gr in strata.items():
    for g, m in gr.items():
        print(f"  [{col} = {g}] studies {m['n_studies']}, lesions {m['n_ref_lesions']}, sensitivity {f(m['lesion_sens'])}, "
              f"FP per scan {f(m['fp_per_scan'])}")
for w in warn: print("  WARNING:", w)
if len(pend): print(f"  WARNING: {len(pend)} object(s) still pending adjudication - left out")
print(f"-> review/evaluation_{tag}.json, _patients.csv, _froc.csv")
```

## tests/check_physics.py

```python
"""v4.18 physics-based synthetic lesions: dipole field (css_common.susceptibility_signal) against analytic values,
then make_synthetic.py --physics on PH1: blooming grows with TE, calcium has the opposite phase, lesions detected.
usage: check_physics.py    (needs PH1 from phantom.py, aligned, in $CSS_BASE)"""
import os, sys, re, shutil, subprocess, numpy as np, nibabel as nib, pandas as pd
B = os.environ["CSS_BASE"]; S = f"{B}/scripts"; py = sys.executable
sys.path.insert(0, S)
from css_common import susceptibility_signal, GAMMA_HZ_T
def fail(m): print("FAIL:", m); sys.exit(1)
def run(*a):
    r = subprocess.run([py, *a], capture_output=True, text=True)
    if r.returncode: print(r.stdout[-600:], r.stderr[-1500:]); fail("step failed: " + " ".join(a))
    return r.stdout

# 1. analytic checks: long cylinder inside field = chi (3 cos^2 th - 1) / 6; sphere outside: + along B0, - across
sh, vox, aff = (80, 80, 80), (0.5, 0.5, 0.5), np.diag([0.5, 0.5, 0.5, 1.0])
k = 2 * np.pi * GAMMA_HZ_T * 3.0 * 0.02 * 1e-6
for name, sl, th in (("parallel", (slice(39, 41), slice(39, 41), slice(5, 75)), 0.0),
                     ("perpendicular", (slice(5, 75), slice(39, 41), slice(39, 41)), np.pi / 2)):
    chi = np.zeros(sh); chi[sl] = 0.45
    ph, _ = susceptibility_signal(chi, vox, aff, 3.0, 0.02)
    want = k * 0.45 * (3 * np.cos(th) ** 2 - 1) / 6; got = float(np.median(ph[sl]))
    if abs(got - want) > 0.1 * abs(want): fail(f"cylinder {name}: phase {got:.2f} vs analytic {want:.2f}")
x, y, z = np.meshgrid(*[np.arange(80)] * 3, indexing="ij")
chi = np.zeros(sh); chi[(x - 40) ** 2 + (y - 40) ** 2 + (z - 40) ** 2 < 9] = 1.0
ph, _ = susceptibility_signal(chi, vox, aff, 3.0, 0.02)
if not (ph[40, 40, 46] > 0 > ph[46, 40, 40] and abs(ph[40, 40, 46] / ph[46, 40, 40] + 2) < 0.4):
    fail(f"sphere dipole pattern {ph[40, 40, 46]:.2f} / {ph[46, 40, 40]:.2f} (expect ratio -2)")
# B0 follows the affine: the same cylinder stored with axes swapped (array axis 0 = scanner z) reads parallel
aff2 = np.zeros((4, 4)); aff2[0, 2] = aff2[1, 1] = aff2[2, 0] = 0.5; aff2[3, 3] = 1
chi = np.zeros(sh); chi[5:75, 39:41, 39:41] = 0.45
ph, _ = susceptibility_signal(chi, vox, aff2, 3.0, 0.02)
if abs(float(np.median(ph[5:75, 39:41, 39:41])) - k * 0.45 / 3) > 0.1 * k * 0.45 / 3: fail("B0 direction not taken from the affine")
print(f"   dipole model: cylinders match the analytic field (within 10 %), sphere pattern +2:-1, B0 from the affine")

# 2. make_synthetic --physics on PH1: TE dependence, calcium phase sign, detection
shutil.copy(f"{B}/synthseg/PH1_seg.nii.gz", f"{B}/synthseg/PH1S_seg.nii.gz")
I = nib.load(f"{B}/data/PH1_swi.nii").get_fdata(); res = {}
for te in (10, 40):
    run(f"{S}/make_synthetic.py", "PH1", "8", "0.6", "1", "--physics", "--te-ms", str(te), "--calcium", "2", "--veins", "4")
    Sx = nib.load(f"{B}/data/PH1S_swi.nii").get_fdata(); t = nib.load(f"{B}/work/PH1S_truth.nii.gz").get_fdata() > 0
    v = nib.load(f"{B}/work/PH1S_veins.nii.gz").get_fdata().astype(int); ph = nib.load(f"{B}/work/PH1S_dphase.nii.gz").get_fdata()
    kd = pd.read_csv(f"{B}/work/PH1S_vein_kinds.csv").set_index("vein_id").kind
    ca = np.isin(v, [i for i, x in kd.items() if x == "calcium"])
    if not ca.any(): fail("no calcium decoys inserted")
    res[te] = dict(bloom=int(((Sx < 0.8 * I) & (I > 0) & ~t & (v == 0)).sum()), ca_phase=float(np.median(ph[ca])),
                   lesion=float(np.median(Sx[t] / I[t])))
if not res[40]["bloom"] > 1.5 * res[10]["bloom"]: fail(f"blooming does not grow with TE: {res}")
if not (res[10]["ca_phase"] < 0 and res[40]["lesion"] < res[10]["lesion"]): fail(f"calcium phase / TE darkening wrong: {res}")
run(f"{S}/make_synthetic.py", "PH1", "8", "0.6", "1", "--physics", "--calcium", "2", "--veins", "4", "--surface-veins", "2")
run(f"{S}/align_seg.py", "PH1S"); run(f"{S}/detect_css.py", "PH1S")
r = re.search(r"found=(\d+) total=(\d+)", run(f"{S}/eval_synthetic.py", "PH1S"))
if int(r.group(1)) < 6: fail(f"physics lesions: only {r.group(1)}/{r.group(2)} found")
print(f"   physics lesions on PH1: darkened voxels outside objects TE 10 ms {res[10]['bloom']} -> TE 40 ms {res[40]['bloom']}; "
      f"calcium phase {res[10]['ca_phase']:+.2f} rad (diamagnetic); found {r.group(1)}/{r.group(2)} at TE 20 ms")
print("RESULT physics OK")
```

## tests/check_evaluate.py

```python
"""v4.18 evaluate.py on hand-made label tables with known answers (no images needed).
usage: check_evaluate.py    (writes into $CSS_BASE/eval_check/)"""
import os, sys, json, subprocess, pandas as pd, numpy as np
B0 = os.environ["CSS_BASE"]; B = f"{B0}/eval_check"; S = f"{B0}/scripts"
def fail(m): print("FAIL:", m); sys.exit(1)
for d in ("review", "labels"): os.makedirs(f"{B}/{d}", exist_ok=True)
def cands(s, n, off=0):
    return pd.DataFrame(dict(cand_id=range(1, n + 1), vox_i=[10 * i + off for i in range(1, n + 1)], vox_j=20, vox_k=5,
                             volume_mm3=30.0, score=np.linspace(2, 1, n)))
uid = lambda s, r: f"{s}:{int(r.vox_i)}:{int(r.vox_j)}:{int(r.vox_k)}:{float(r.volume_mm3):g}"
cons = []
def lab(s, r, cls, lst="candidate", oid=None):
    cons.append(dict(study_id=s, cand_uid=uid(s, r), list=lst, obj_id=oid or int(r.get("cand_id", 0) or 0), readers="A+B",
                     adjudicator="", consensus_label=cls, consensus_class=cls, method="agreement"))
# S1 (test, positive): ranks 1 and 4 cSS, 2 vein, 3 uncertain, 5 unlabelled; one EXCLUDED object called cSS
c1 = cands("S1", 5); c1.to_csv(f"{B}/review/S1_candidates.csv", index=False)
x1 = cands("S1", 1, off=3).rename(columns={"cand_id": "excl_id"}).assign(rule_version="E7@v4.12")
x1.to_csv(f"{B}/review/S1_excluded.csv", index=False)
for rk, cls in ((1, "cSS"), (2, "other_vein"), (3, "uncertain"), (4, "cSS")): lab("S1", c1.iloc[rk - 1], cls)
lab("S1", x1.iloc[0], "cSS", "excluded", 1)
# S2 (test, negative): 3 candidates, all veins; plus a stale label (no matching object)
c2 = cands("S2", 3); c2.to_csv(f"{B}/review/S2_candidates.csv", index=False)
pd.DataFrame(columns=["excl_id"]).to_csv(f"{B}/review/S2_excluded.csv", index=False)
for rk in (1, 2, 3): lab("S2", c2.iloc[rk - 1], "other_vein")
cons.append(dict(study_id="S2", cand_uid="S2:1:1:1:9", list="candidate", obj_id=9, readers="A+B", adjudicator="",
                 consensus_label="cSS", consensus_class="cSS", method="agreement"))
# S3 (dev): must not appear in the test evaluation
c3 = cands("S3", 2); c3.to_csv(f"{B}/review/S3_candidates.csv", index=False)
lab("S3", c3.iloc[0], "cSS")
pd.DataFrame(cons).to_csv(f"{B}/labels/consensus.csv", index=False)
# missed outlines in S1: one never generated (two readers, 2 voxels apart = one lesion), one on shown rank 2 (the 'vein')
pd.DataFrame([dict(study_id="S1", reader_id="A", vox_i=100, vox_j=100, vox_k=5, fate="not generated", obj_id=np.nan, rule_version="", reason="not dark"),
              dict(study_id="S1", reader_id="B", vox_i=101, vox_j=101, vox_k=5, fate="not generated", obj_id=np.nan, rule_version="", reason="not dark"),
              dict(study_id="S1", reader_id="A", vox_i=20, vox_j=20, vox_k=5, fate="shown", obj_id=2, rule_version="", reason="")]
             ).to_csv(f"{B}/labels/missed_lesions.csv", index=False)
pd.DataFrame([dict(study_id="S1", css_status="pos", split="test", slice_mm="1.2"),
              dict(study_id="S2", css_status="neg", split="test", slice_mm="2.0"),
              dict(study_id="S3", css_status="pos", split="dev", slice_mm="1.2")]).to_csv(f"{B}/labels/study_meta.csv", index=False)
r = subprocess.run([sys.executable, f"{S}/evaluate.py", "--k", "1,2", "--boot", "200", "--by", "slice_group"],
                   capture_output=True, text=True, env=dict(os.environ, CSS_BASE=B))
if r.returncode: print(r.stdout, r.stderr); fail("evaluate.py failed")
o = json.load(open(f"{B}/review/evaluation_test_consensus.json"))["overall"]
want = dict(n_studies=2, n_ref_lesions=5, lesion_sens=0.6, sens_top1=0.2, sens_top2=0.4, patient_sens=1.0,
            fp_per_scan=1.5, fp_per_neg_scan=3.0, uncertain_per_scan=0.5, fp_per_scan_top2=1.0)
bad = {k: (o[k], v) for k, v in want.items() if abs(o[k] - v) > 1e-9}
if bad: print(r.stdout); fail(f"metrics (got, want): {bad}")
if "S2: 1 label(s) match no current" not in r.stdout or "S1: 1 shown candidate(s) without a label" not in r.stdout:
    fail("stale / unlabelled warnings missing")
fr = pd.read_csv(f"{B}/review/evaluation_test_consensus_froc.csv")
if not (fr.lesion_sens.is_monotonic_increasing and fr.fp_per_scan.is_monotonic_increasing): fail("FROC not monotone")
print("   " + "\n   ".join(l for l in r.stdout.splitlines()[:6]))
print("RESULT evaluate OK")
```

## scripts/mp2rage_robust.py

```python
#!/usr/bin/env python3
"""MP2RAGE background clean-up before FreeSurfer (v4.18): the 'robust combination' of O'Brien et al. 2014
(PLoS One 9:e99676) as implemented by J. Marques (RobustCombination.m), from MAGNITUDE images.
usage: mp2rage_robust.py UNI.nii.gz INV1.nii.gz INV2.nii.gz OUT_T1.nii.gz [--lambda auto|1|5|...]
The MP2RAGE UNI image has salt-and-pepper noise outside the head (division of two noisy signals). recon-all
then fails to strip dura / background (LEMON P201: 'extra large defect ... dura has not been removed from
wm.mgz', topology fixing aborted). The robust combination adds a noise term beta to the UNI formula:
      UNI_r = (INV1 * INV2 - beta) / (INV1^2 + INV2^2 + 2 beta),   beta = lambda * mean(INV2 noise)^2
which leaves the brain (high INV2 signal) unchanged and pushes the background to -0.5 (black). INV1's sign
(lost in magnitude images) is recovered from UNI: INV1 = the root of UNI = INV1 INV2 / (INV1^2 + INV2^2) that
is closer to the measured |INV1|. --lambda auto (v4.18.2): the smallest of 1, 2, 4 ... 128 that leaves < 5 % of the
background band 5-20 mm outside the head brighter than 25 % of the brain, without changing the brain (r >= 0.95). Output is scaled to 0..4095 like a scanner UNI. Prints the
correlation with the original UNI inside the head (should stay > 0.95)."""
import sys, argparse, numpy as np, nibabel as nib
from scipy import ndimage as ndi

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("uni"); ap.add_argument("inv1"); ap.add_argument("inv2"); ap.add_argument("out")
ap.add_argument("--lambda", dest="lam", default="auto")
a = ap.parse_args()
iu, i1, i2 = (nib.load(p) for p in (a.uni, a.inv1, a.inv2))
if not (iu.shape[:3] == i1.shape[:3] == i2.shape[:3]) or not (np.allclose(iu.affine, i1.affine, atol=1e-2)
                                                             and np.allclose(iu.affine, i2.affine, atol=1e-2)):
    sys.exit("UNI, INV1 and INV2 must be on the same grid (same MP2RAGE acquisition)")
U = np.asarray(iu.dataobj, np.float64); A1 = np.abs(np.asarray(i1.dataobj, np.float64)); A2 = np.abs(np.asarray(i2.dataobj, np.float64))
# UNI to [-0.5, 0.5]
if U.max() > 1.5:
    lo, hi = (0.0, 4095.0) if U.min() >= 0 and U.max() <= 4096 else (float(U.min()), float(U.max()))
    u = (U - lo) / (hi - lo) - 0.5
else:
    u = U - 0.5 if U.min() >= 0 else U
u = np.clip(u, -0.5, 0.5)
# INV1 sign from UNI: u (i1^2 + i2^2) = i1 i2  ->  u i1^2 - i2 i1 + u i2^2 = 0
with np.errstate(invalid="ignore", divide="ignore"):
    disc = np.sqrt(np.clip(A2 ** 2 - 4 * u * u * A2 ** 2, 0, None))
    r_pos = np.where(np.abs(u) > 1e-6, (A2 + disc) / (2 * u), 0.0)
    r_neg = np.where(np.abs(u) > 1e-6, (A2 - disc) / (2 * u), 0.0)
I1 = np.where(np.abs(np.abs(r_pos) - A1) < np.abs(np.abs(r_neg) - A1), r_pos, r_neg)
I1 = np.where(np.isfinite(I1), I1, 0.0)
# v4.18.2: noise from the background band 5-20 mm OUTSIDE the head (INV2 has a clean head outline). The image border
# is often zero-padded / defaced (LEMON: border noise 10 -> lambda 1, speckle left around the skull and neck).
from skimage.filters import threshold_otsu
vox = np.array(iu.header.get_zooms()[:3], float)
pos = A2[A2 > 0]
head = A2 > threshold_otsu(pos) * 0.5
head = ndi.binary_closing(head, iterations=2); head = ndi.binary_fill_holes(head)
lab, n = ndi.label(head)
if n: head = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
d_out = ndi.distance_transform_edt(~head, sampling=vox)
shell = (d_out > 5) & (d_out <= 20) & (A2 > 0)
if shell.sum() < 1000:                                    # tiny field of view: fall back to everything outside
    shell = (d_out > 2) & (A2 > 0)
noise = float(np.mean(A2[shell])) if shell.any() else float(np.percentile(pos, 1))
core = ndi.binary_erosion(head, iterations=max(1, int(round(10 / vox.min()))))
if not core.any(): core = head
ref = (u[core] + 0.5) * 4095
def robust(lam):
    b = lam * noise ** 2
    return (I1 * A2 - b) / (I1 ** 2 + A2 ** 2 + 2 * b)
lams = [1, 2, 4, 8, 16, 32, 64, 128] if a.lam == "auto" else [float(a.lam)]
prev = None
for lam in lams:
    out_l = np.clip((robust(lam) + 0.5) * 4095, 0, 4095)
    hm = float(np.median(out_l[core]))
    speck = float((out_l[shell] > 0.25 * hm).mean()) if shell.any() else 0.0
    rc = float(np.corrcoef(out_l[core], ref)[0, 1])
    if rc < 0.95 and prev is not None:                    # brain contrast would change: keep the previous lambda
        lam, out, speck, ok = prev; break
    out, ok = out_l, speck < 0.05
    prev = (lam, out_l, speck, ok)
    if ok: break
bd = shell
c = float(np.corrcoef(out[core], ref)[0, 1])
img = nib.Nifti1Image(out.astype(np.float32), iu.affine, iu.header); img.set_data_dtype(np.float32)
nib.save(img, a.out)
print(f"{a.out}: robust MP2RAGE combination, lambda {lam:g} (background noise {noise:.1f}, band 5-20 mm outside the head); "
      f"speckle left in that band {speck:.1%} (was {float(((u[shell] + 0.5) * 4095 > 0.25 * np.median(ref)).mean()):.1%}); "
      f"correlation with the original UNI inside the brain {c:.3f}" + ("" if ok else "  WARNING: background still speckled"))
if c < 0.95: print("WARNING: brain contrast changed (correlation < 0.95) - check the image before recon-all")
```

## tests/check_mp2rage.py

```python
"""v4.18 mp2rage_robust.py: synthetic MP2RAGE (signed INV1, complex noise) -> background noise removed, brain kept.
usage: check_mp2rage.py    (scratch files in $CSS_BASE/work)"""
import os, sys, subprocess, numpy as np, nibabel as nib
B = os.environ["CSS_BASE"]; W = f"{B}/work"; os.makedirs(W, exist_ok=True)
def fail(m): print("FAIL:", m); sys.exit(1)
rng = np.random.default_rng(0); sh = (100, 110, 90)
x, y, z = np.meshgrid(*[np.arange(n) for n in sh], indexing="ij")
e = lambda a, b, c: ((x - 50) / a) ** 2 + ((y - 55) / b) ** 2 + ((z - 45) / c) ** 2 < 1
head, brain, wm = e(40, 45, 38), e(32, 37, 30), e(20, 24, 18)
i1 = np.where(head, -200.0, 0.0); i1[brain] = -80; i1[wm] = 150          # INV1 sign flips with T1
i2 = np.where(head, 600.0, 0.0); i2[brain] = 700; i2[wm] = 800
n = lambda v: v + rng.normal(0, 20, sh) + 1j * rng.normal(0, 20, sh)
c1, c2 = n(i1), n(i2)
uni = np.real(c1 * np.conj(c2)) / (np.abs(c1) ** 2 + np.abs(c2) ** 2)
pad = np.zeros(sh, bool); pad[8:-8, 8:-8, 8:-8] = True       # LEMON-like zero-padded / defaced border (v4.18.2)
for name, v in (("uni", np.round((uni + 0.5) * 4095)), ("inv1", np.abs(c1)), ("inv2", np.abs(c2))):
    nib.save(nib.Nifti1Image((v * pad).astype(np.float32), np.eye(4)), f"{W}/MP2_{name}.nii.gz")
r = subprocess.run([sys.executable, f"{B}/scripts/mp2rage_robust.py", *(f"{W}/MP2_{k}.nii.gz" for k in ("uni", "inv1", "inv2")),
                    f"{W}/MP2_robust.nii.gz"], capture_output=True, text=True)
if r.returncode: print(r.stdout, r.stderr); fail("mp2rage_robust failed")
U = nib.load(f"{W}/MP2_uni.nii.gz").get_fdata(); O = nib.load(f"{W}/MP2_robust.nii.gz").get_fdata()
from scipy import ndimage as ndi
band = (ndi.distance_transform_edt(~head) > 5) & (ndi.distance_transform_edt(~head) <= 15) & pad
hm = np.median(U[wm])
bg_before, bg_after = (U[band] > 0.25 * hm).mean(), (O[band] > 0.25 * hm).mean()
cor = float(np.corrcoef(U[brain], O[brain])[0, 1])
if not (bg_after < 0.05 and bg_before > 0.3 and cor > 0.95 and abs(np.median(O[wm]) - np.median(U[wm])) < 0.10 * np.median(U[wm])):
    fail(f"speckle {bg_before:.2f} -> {bg_after:.3f}, brain r {cor:.3f}, wm {np.median(O[wm]):.0f} vs {np.median(U[wm]):.0f}; {r.stdout.strip()}")
print(f"   robust MP2RAGE (zero-padded border): speckle 5-15 mm outside the head {bg_before:.0%} -> {bg_after:.1%}, "
      f"brain correlation {cor:.3f}, white matter within 10 %; {r.stdout.strip().split('; ')[0].split(': ')[-1]}")
print("RESULT mp2rage OK")
```
