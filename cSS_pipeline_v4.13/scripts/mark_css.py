#!/usr/bin/env python3
"""Mark reviewed candidates and score.
usage: mark_css.py SUBJECT 1,4,9 [--ich 3,7] [--reviewed 1-40 list]     (use 'none' when no cSS)
--reviewed: ids the reader actually looked at (review_css.py passes this). Default: all candidates.
--ich: candidates you judged to be siderosis/hemosiderin CONNECTED TO a lobar ICH.
       They are excluded from the multifocality score (standard convention) but reported."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, subprocess, pandas as pd
if len(sys.argv) < 3:
    sys.exit("usage: mark_css.py SUBJECT 1,4,9 [--ich 3,7]   (or: none)")
subj, arg = sys.argv[1], sys.argv[2]
ich_arg = sys.argv[sys.argv.index("--ich") + 1] if "--ich" in sys.argv else ""
rev_arg = sys.argv[sys.argv.index("--reviewed") + 1] if "--reviewed" in sys.argv else None
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
csv = f"{base}/review/{subj}_candidates.csv"
df = pd.read_csv(csv)
parse = lambda a: [] if a.lower() in ("", "none") else [int(x) for x in a.split(",") if x.strip()]
ids, ich = parse(arg), parse(ich_arg)
rev = parse(rev_arg) if rev_arg is not None else df.cand_id.astype(int).tolist()
valid = set(df.cand_id)
bad = [i for i in ids + ich + rev if i not in valid]
if bad: sys.exit(f"candidate ids not in list: {bad} (valid 1-{df.cand_id.max()})")
both = sorted(set(ids) & set(ich))
if both: sys.exit(f"ids marked both cSS and near-ICH: {both}")
df["accept"] = df.cand_id.isin(ids).astype(int)
df["near_ich"] = df.cand_id.isin(ich).astype(int)
df["reviewed"] = df.cand_id.isin(set(rev) | set(ids) | set(ich)).astype(int)
df.to_csv(csv, index=False)
print(f"{subj}: accepted {len(ids)} of {len(df)} candidates {ids if ids else ''}"
      + (f"; near-ICH (not scored) {ich}" if ich else ""))
subprocess.run([sys.executable, f"{base}/scripts/score_css.py", subj], check=True)
