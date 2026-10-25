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
    # v5.1: what the detector did with every generated object in this study (review/ID_manifest.json fates)
    _mf = f"{B}/review/{s}_manifest.json"
    _f = (json.load(open(_mf)).get("detect", {}).get("fates") or {}) if os.path.exists(_mf) else {}
    row.update(obj_generated=_f.get("objects"), obj_dropped=sum((_f.get("dropped_by_gate") or {}).values()) if _f else None,
               obj_excluded=sum((_f.get("excluded_by_rule") or {}).values()) if _f else None,
               obj_shown=_f.get("shown"), obj_qc_blocked=_f.get("qc_blocked"))
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
