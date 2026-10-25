"""v4.18 evaluate.py on hand-made label tables with known answers (no images needed).
usage: check_evaluate.py    (writes into $CSS_BASE/eval_check/)"""
import os, sys, json, subprocess, pandas as pd, numpy as np
B0 = os.environ["CSS_BASE"]; B = f"{B0}/eval_check"; S = f"{B0}/scripts"
def fail(m): print("FAIL:", m); sys.exit(1)
for d in ("review", "labels"): os.makedirs(f"{B}/{d}", exist_ok=True)
def cands(s, n, off=0):
    return pd.DataFrame(dict(cand_id=range(1, n + 1), vox_i=[10 * i + off for i in range(1, n + 1)], vox_j=20, vox_k=5,
                             volume_mm3=30.0, score=np.linspace(2, 1, n)))
uid = lambda s, r: f"{s}:{int(r.vox_i)}:{int(r.vox_j)}:{int(r.vox_k)}:{float(r.volume_mm3):g}"
cons = []
def lab(s, r, cls, lst="candidate", oid=None):
    cons.append(dict(study_id=s, cand_uid=uid(s, r), list=lst, obj_id=oid or int(r.get("cand_id", 0) or 0), readers="A+B",
                     adjudicator="", consensus_label=cls, consensus_class=cls, method="agreement"))
# S1 (test, positive): ranks 1 and 4 cSS, 2 vein, 3 uncertain, 5 unlabelled; one EXCLUDED object called cSS
c1 = cands("S1", 5); c1.to_csv(f"{B}/review/S1_candidates.csv", index=False)
x1 = cands("S1", 1, off=3).rename(columns={"cand_id": "excl_id"}).assign(rule_version="E7@v4.12")
x1.to_csv(f"{B}/review/S1_excluded.csv", index=False)
for rk, cls in ((1, "cSS"), (2, "other_vein"), (3, "uncertain"), (4, "cSS")): lab("S1", c1.iloc[rk - 1], cls)
lab("S1", x1.iloc[0], "cSS", "excluded", 1)
# S2 (test, negative): 3 candidates, all veins; plus a stale label (no matching object)
c2 = cands("S2", 3); c2.to_csv(f"{B}/review/S2_candidates.csv", index=False)
pd.DataFrame(columns=["excl_id"]).to_csv(f"{B}/review/S2_excluded.csv", index=False)
for rk in (1, 2, 3): lab("S2", c2.iloc[rk - 1], "other_vein")
cons.append(dict(study_id="S2", cand_uid="S2:1:1:1:9", list="candidate", obj_id=9, readers="A+B", adjudicator="",
                 consensus_label="cSS", consensus_class="cSS", method="agreement"))
# S3 (dev): must not appear in the test evaluation
c3 = cands("S3", 2); c3.to_csv(f"{B}/review/S3_candidates.csv", index=False)
lab("S3", c3.iloc[0], "cSS")
pd.DataFrame(cons).to_csv(f"{B}/labels/consensus.csv", index=False)
# missed outlines in S1: one never generated (two readers, 2 voxels apart = one lesion), one on shown rank 2 (the 'vein')
pd.DataFrame([dict(study_id="S1", reader_id="A", vox_i=100, vox_j=100, vox_k=5, fate="not generated", obj_id=np.nan, rule_version="", reason="not dark"),
              dict(study_id="S1", reader_id="B", vox_i=101, vox_j=101, vox_k=5, fate="not generated", obj_id=np.nan, rule_version="", reason="not dark"),
              dict(study_id="S1", reader_id="A", vox_i=20, vox_j=20, vox_k=5, fate="shown", obj_id=2, rule_version="", reason="")]
             ).to_csv(f"{B}/labels/missed_lesions.csv", index=False)
pd.DataFrame([dict(study_id="S1", css_status="pos", split="test", slice_mm="1.2"),
              dict(study_id="S2", css_status="neg", split="test", slice_mm="2.0"),
              dict(study_id="S3", css_status="pos", split="dev", slice_mm="1.2")]).to_csv(f"{B}/labels/study_meta.csv", index=False)
r = subprocess.run([sys.executable, f"{S}/evaluate.py", "--k", "1,2", "--boot", "200", "--by", "slice_group"],
                   capture_output=True, text=True, env=dict(os.environ, CSS_BASE=B))
if r.returncode: print(r.stdout, r.stderr); fail("evaluate.py failed")
o = json.load(open(f"{B}/review/evaluation_test_consensus.json"))["overall"]
want = dict(n_studies=2, n_ref_lesions=5, lesion_sens=0.6, sens_top1=0.2, sens_top2=0.4, patient_sens=1.0,
            fp_per_scan=1.5, fp_per_neg_scan=3.0, uncertain_per_scan=0.5, fp_per_scan_top2=1.0)
bad = {k: (o[k], v) for k, v in want.items() if abs(o[k] - v) > 1e-9}
if bad: print(r.stdout); fail(f"metrics (got, want): {bad}")
if "S2: 1 label(s) match no current" not in r.stdout or "S1: 1 shown candidate(s) without a label" not in r.stdout:
    fail("stale / unlabelled warnings missing")
fr = pd.read_csv(f"{B}/review/evaluation_test_consensus_froc.csv")
if not (fr.lesion_sens.is_monotonic_increasing and fr.fp_per_scan.is_monotonic_increasing): fail("FROC not monotone")
print("   " + "\n   ".join(l for l in r.stdout.splitlines()[:6]))
print("RESULT evaluate OK")
