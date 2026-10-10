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
