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
