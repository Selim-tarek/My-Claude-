#!/usr/bin/env python3
"""Import an expert's answers from the blinded sheet (expert_sheet.py) as REFERENCE calls.
usage: expert_import.py SUBJECT READER --css 3,17,40 [--vein 1,2] [--ich 9] [--unsure 5]
       expert_import.py SUBJECT READER --letters "1V 2V 3C 4N 5U ..."      (letters from the grid)
Sheet numbers not listed are recorded as N (normal / other).
Writes review/ID_expert_READER.csv (cand_id, call, fingerprint) - the reader's own calls in
review/ID_calls.csv are left untouched. Then scores the expert calls:
  expert_import.py ... --score     (runs mark_css.py + score_css.py with the expert's cSS / ICH calls)"""
import sys, os, re, subprocess, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, FP, match_calls

if len(sys.argv) < 3: sys.exit(__doc__)
subj, reader = sys.argv[1], sys.argv[2]
base = _base()
key = pd.read_csv(f"{base}/review/{subj}_expert_key.csv")
arg = lambda f: sys.argv[sys.argv.index(f) + 1] if f in sys.argv else ""
nums = lambda s: [int(x) for x in re.split(r"[,\s]+", s.strip()) if x]
NAME = {"C": "cSS", "V": "Vein", "N": "Normal", "H": "Near ICH", "U": "Unsure"}
calls = {n: "Normal" for n in key.sheet_no}
if arg("--letters"):
    for tok in arg("--letters").split():
        m = re.fullmatch(r"(\d+)([CVNHU])", tok.strip().upper())
        if not m: sys.exit(f"cannot read '{tok}' (expected e.g. 12C)")
        calls[int(m.group(1))] = NAME[m.group(2)]
for flag, call in (("--css", "cSS"), ("--vein", "Vein"), ("--ich", "Near ICH"), ("--unsure", "Unsure")):
    for n in nums(arg(flag)): calls[n] = call
bad = [n for n in calls if n not in set(key.sheet_no)]
if bad: sys.exit(f"sheet numbers not on the sheet: {bad}")
key["call"] = key.sheet_no.map(calls)
# re-attach to the CURRENT candidates by fingerprint (the detector may have been re-run)
cand = pd.read_csv(f"{base}/review/{subj}_candidates.csv")
mapped, dropped = match_calls(cand, key[["cand_id", "call"] + FP])
if dropped: print(f"WARNING: {dropped} answers match no current candidate (detector re-run since the sheet)")
c = cand.set_index("cand_id")
rows = [dict(cand_id=k, call=v, **{f: c.loc[k, f] for f in FP}) for k, v in sorted(mapped.items())]
out = f"{base}/review/{subj}_expert_{reader}.csv"
pd.DataFrame(rows, columns=["cand_id", "call"] + FP).to_csv(out, index=False)
summ = pd.Series(list(mapped.values())).value_counts().to_dict()
print(f"{subj}: expert {reader} -> {out}   {summ}")
if "--score" in sys.argv:
    css = sorted(k for k, v in mapped.items() if v == "cSS"); ich = sorted(k for k, v in mapped.items() if v == "Near ICH")
    cmd = [sys.executable, f"{base}/scripts/mark_css.py", subj, ",".join(map(str, css)) or "none"]
    if ich: cmd += ["--ich", ",".join(map(str, ich))]
    subprocess.run(cmd, check=True)
