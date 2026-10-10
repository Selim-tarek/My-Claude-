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
  not generated  never a dark ridge component (contrast / ridge threshold / outside the search zone)
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
        r["fate"] = "not generated"
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
