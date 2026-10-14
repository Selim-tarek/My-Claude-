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
