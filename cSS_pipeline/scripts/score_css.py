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
