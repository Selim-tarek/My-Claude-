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
