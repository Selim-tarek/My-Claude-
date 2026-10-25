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
