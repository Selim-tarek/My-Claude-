#!/usr/bin/env python3
"""Per-case audit report (v5.1): what went in, what the detector did with every object, and what to distrust.
usage: audit_case.py SUBJECT            -> prints the report, writes review/ID_audit.json
Reads only existing outputs (review/ID_manifest.json, ID_qc.json, candidates / excluded CSVs, work/ID_dropped.csv)
and image headers; changes nothing. Sections:
  inputs       file, sha256 (from the manifest), shape, voxel size, orientation, obliquity, S-I array axis
  acquisition  sequence type (SWI / T2*-GRE / SWI built from GRE by make_swi.py / unknown), field, TE, TR, slice
               thickness, vendor (dcm2niix JSON, if present)
  anatomy      T1 route, FLAIR / phase availability and FRESHNESS against the SWI, registration method, cost,
               label shift, edge agreement, QC status and reasons
  fates        raw components -> not generated (no strong seed) -> pieces -> merged objects -> dropped per gate /
               excluded per rule / shown; the counts must add up (the detector stops otherwise, exit 4)
Research software: expert review of every candidate remains required; not validated for clinical use."""
import sys, os, json, glob, nibabel as nib, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, stale_against

if len(sys.argv) < 2: sys.exit(__doc__)
s = sys.argv[1]; B = _base()
J = lambda p: json.load(open(p)) if os.path.exists(p) else {}
man, qc = J(f"{B}/review/{s}_manifest.json"), J(f"{B}/review/{s}_qc.json")
det = man.get("detect", {})
if not det: sys.exit(f"{s}: no detector run recorded (review/{s}_manifest.json) - run detect_css.py first")
swi_p = f"{B}/data/{s}_swi.nii"
rep = dict(subject=s, schema_version=man.get("schema_version"), build=(det.get("environment") or {}).get("build"),
           detector=det.get("pipeline_version"), run_at=det.get("run_at"), parameters=det.get("parameters"),
           note="research software; expert review of every candidate remains required; not validated for clinical use")
# inputs
def geom(p):
    try: im = nib.load(p)
    except Exception: return None
    a = im.affine; vox = [round(float(v), 3) for v in im.header.get_zooms()[:3]]
    si = int(np.argmax(np.abs(a[2, :3]))); n = a[:3, si] / np.linalg.norm(a[:3, si])
    return dict(shape=list(im.shape[:3]), voxel_mm=vox, orientation="".join(nib.aff2axcodes(a)), si_axis=si,
                obliquity_deg=round(float(np.degrees(np.arccos(min(1.0, abs(n[2]))))), 1))
inp = det.get("inputs", {})
rep["inputs"] = {k: dict(path=v.get("path"), sha256=(v.get("sha256") or "")[:16], **(geom(v["path"]) or {}))
                 for k, v in inp.items() if k in ("swi", "labels", "t1_labels", "flair", "phase", "destrieux")}
# acquisition
acq = qc.get("acquisition") or {}
if acq.get("sequence"): seq = acq["sequence"]
elif os.path.exists(f"{B}/work/{s}_gre_swi.nii.gz"): seq = "SWI built from GRE magnitude + phase (make_swi.py)"
else: seq = "unknown (no dcm2niix JSON in raw/ID/SWI)"
rep["acquisition"] = dict(sequence=seq, field_T=acq.get("field_T"), TE_ms=acq.get("TE_ms"), TR_ms=acq.get("TR_ms"),
                          slice_thickness_mm=acq.get("slice_thickness_mm"), vendor=acq.get("manufacturer"))
# anatomy, freshness
fresh = {}
for name, p in (("T1 labels", f"{B}/work/{s}_t1seg_swispace.nii.gz"), ("FLAIR", f"{B}/work/{s}_flair_swispace.nii.gz"),
                ("phase", f"{B}/data/{s}_phase.nii.gz"), ("Destrieux labels", f"{B}/work/{s}_a2009s_swispace.nii.gz")):
    fresh[name] = "absent" if not os.path.exists(p) else ("STALE (older than the SWI) - ignored" if stale_against(p, [swi_p]) else "fresh")
sp = f"{B}/review/{s}_surface_features.csv"
fresh["surface features"] = "absent" if not os.path.exists(sp) else \
    ("STALE (older than this detection) - ignored" if stale_against(sp, [f"{B}/review/{s}_candidates.nii.gz"]) else "fresh")
reg = qc.get("registration") or {}
rep["anatomy"] = dict(route=qc.get("anatomy"), inputs=fresh, registration=reg.get("method"), registration_cost=reg.get("cost"),
                      label_shift_mm=qc.get("label_shift_mm"), edge_bias_mm=qc.get("edge_bias_mm"),
                      edge_scatter_mm=qc.get("edge_scatter_mm"), edge_n=qc.get("edge_n"),
                      qc_status=qc.get("status"), qc_reasons=qc.get("reasons"), qc_notes=qc.get("notes"),
                      phase_rule=(det.get("parameters") or {}).get("phase"))
# fates
F = det.get("fates") or {}
dp, xp, cp = f"{B}/work/{s}_dropped.csv", f"{B}/review/{s}_excluded.csv", f"{B}/review/{s}_candidates.csv"
nd = len(pd.read_csv(dp)) if os.path.exists(dp) else None
nx = len(pd.read_csv(xp)) if os.path.exists(xp) else None
nc = len(pd.read_csv(cp)) if os.path.exists(cp) else None
check = None
if F.get("objects") is not None and None not in (nd, nx, nc):
    check = "ok" if (F.get("qc_blocked") or nd + nx + nc) == F["objects"] else \
        f"MISMATCH: {F['objects']} objects vs {nd} dropped + {nx} excluded + {nc} shown (files changed after the run?)"
rep["fates"] = dict(F, files=dict(dropped=nd, excluded=nx, shown=nc), consistency=check)
json.dump(rep, open(f"{B}/review/{s}_audit.json", "w"), indent=2, default=str)

# printed report
print(f"{s}: audit (schema {rep['schema_version']}, build {rep['build']}, detector {rep['detector']}, run {rep['run_at']})")
print("  research software - expert review of every candidate remains required; not validated for clinical use")
for k, v in rep["inputs"].items():
    print(f"  input {k:10s} {v.get('orientation', '?')} {v.get('voxel_mm')} mm, obliquity {v.get('obliquity_deg')} deg, "
          f"S-I axis {v.get('si_axis')}, sha256 {v.get('sha256')}")
a = rep["acquisition"]
print(f"  acquisition: {a['sequence']}; field {a['field_T']} T, TE {a['TE_ms']} ms, slice {a['slice_thickness_mm']} mm, vendor {a['vendor']}")
an = rep["anatomy"]
print(f"  anatomy: {an['route']}; " + ", ".join(f"{k} {v}" for k, v in an["inputs"].items()))
print(f"  registration: {an['registration']} (cost {an['registration_cost']}); label shift {an['label_shift_mm']} mm, "
      f"edge bias {an['edge_bias_mm']} mm, scatter {an['edge_scatter_mm']} mm (n {an['edge_n']})")
print(f"  QC: {an['qc_status']}" + (f" - {'; '.join(an['qc_reasons'])}" if an["qc_reasons"] else "")
      + (f" [notes: {'; '.join(an['qc_notes'])}]" if an.get("qc_notes") else ""))
print(f"  phase rule: {an['phase_rule']}")
if F:
    print(f"  fates: {F.get('raw_components')} raw components -> {F.get('not_generated_unseeded')} not generated (no strong seed) "
          f"-> {F.get('pieces_before_merge')} pieces -> {F.get('objects')} objects")
    print(f"         dropped by gate {F.get('dropped_by_gate')}; excluded by rule {F.get('excluded_by_rule')}; "
          f"shown {F.get('shown')}; QC-blocked {F.get('qc_blocked')}  [{check}]")
else:
    print("  fates: not recorded (detector run before v5.1)")
