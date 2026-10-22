"""v5.0 checks: train_ranker.py / rank_candidates.py on synthetic label tables with known answers (no images).
usage: check_ranker.py    (writes into $CSS_BASE/ranker_check/)
1. one feature (darkness_z) separates cSS perfectly -> LOPO AUC 1.0, its weight the largest and positive;
   test-split studies are never trained on; fewer than 3 labelled patients -> refused.
2. labels unrelated to every feature -> LOPO AUC near 0.5 with a bootstrap interval that includes 0.5 (honest).
3. rank_candidates.py: score_model / rank_model appended at the END, every other column unchanged; a training
   study is refused; a feature missing from the study is filled with the training median and named."""
import os, sys, json, subprocess, numpy as np, pandas as pd
B0 = os.environ.get("CSS_BASE", "/tmp"); B = f"{B0}/ranker_check"; S = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts")
py = sys.executable
def fail(m): print("FAIL:", m); sys.exit(1)
def run(*a, ok=True, env=None):
    e = dict(os.environ, CSS_BASE=B); e.update(env or {})
    r = subprocess.run([py, *a], capture_output=True, text=True, env=e)
    if ok and r.returncode: print(r.stdout[-1500:], r.stderr[-1500:]); fail("step failed: " + " ".join(a))
    return r
sys.path.insert(0, S)
from css_common import cand_uid

def make(signal, seed, n_pat=6, n_cand=30, test_pat=1):
    rng = np.random.default_rng(seed)
    for d in ("review", "labels"): os.makedirs(f"{B}/{d}", exist_ok=True)
    for f in os.listdir(f"{B}/review"): os.remove(f"{B}/review/{f}")
    meta, lab = [], []
    for p in range(n_pat + test_pat):
        s = f"Q{p:02d}"; y = (rng.random(n_cand) < 0.3).astype(int); y[0], y[1] = 1, 0
        dz = (2.5 + 2.0 * y + rng.normal(0, 0.3, n_cand)) if signal else rng.normal(3, 1, n_cand)
        c = pd.DataFrame(dict(cand_id=np.arange(1, n_cand + 1), vox_i=rng.integers(0, 200, n_cand), vox_j=rng.integers(0, 200, n_cand),
                              vox_k=np.arange(n_cand), volume_mm3=np.round(rng.uniform(25, 400, n_cand), 1),
                              darkness_z=np.round(dz, 3), tube_ratio=np.round(rng.uniform(0, 1, n_cand), 3),
                              region_z=np.round(rng.normal(-3, 1, n_cand), 3), score_v4=np.round(rng.uniform(0, 3, n_cand), 3),
                              vein_tree_mm=255.0, mirror_dark_frac=rng.uniform(0, 1, n_cand)))
        c.to_csv(f"{B}/review/{s}_candidates.csv", index=False)
        meta.append(dict(study_id=s, split="test" if p >= n_pat else "dev", css_status="pos"))
        for (_, r), yy in zip(c.iterrows(), y):
            lab.append(dict(study_id=s, cand_uid=cand_uid(s, r), consensus_class="cSS" if yy else "surface_vein", method="agree"))
    pd.DataFrame(meta).to_csv(f"{B}/labels/study_meta.csv", index=False)
    pd.DataFrame(lab).to_csv(f"{B}/labels/consensus.csv", index=False)

# 1. perfect signal
make(True, 1)
out = run(f"{S}/train_ranker.py", "--boot", "200").stdout
M = json.load(open(f"{B}/review/ranker_model.json"))
w = dict(zip(M["features"], M["weights"]))
if M["lopo_auc"] != 1.0: fail(f"perfect feature: LOPO AUC {M['lopo_auc']}")
if max(w, key=lambda k: abs(w[k])) != "darkness_z" or w["darkness_z"] <= 0: fail(f"weights {w}")
if any(s.startswith("Q06") for s in M["training_studies"]): fail("trained on a test-split study")
if "vein_tree_mm" in M["features"] or "mirror_dark_frac" in M["features"]: fail("deprecated column used as a feature")
print(f"perfect feature: LOPO AUC {M['lopo_auc']} {M['lopo_auc_ci95']}, darkness_z weight {w['darkness_z']:+.2f} (largest); "
      f"score_v4 baseline {M['baseline_auc']}; trained on {len(M['training_studies'])} dev studies, test study left out")
r = run(f"{S}/train_ranker.py", "--split", "test", ok=False)
if r.returncode == 0: fail("--split test was accepted")
# 3. rank_candidates
before = pd.read_csv(f"{B}/review/Q06_candidates.csv")
pd.read_csv(f"{B}/review/Q06_candidates.csv").drop(columns="region_z").to_csv(f"{B}/review/Q06_candidates.csv", index=False)
before = before.drop(columns="region_z")
o = run(f"{S}/rank_candidates.py", "Q06").stdout
after = pd.read_csv(f"{B}/review/Q06_candidates.csv")
if list(after.columns) != list(before.columns) + ["score_model", "rank_model"]: fail(f"columns {list(after.columns)}")
if not after[before.columns].equals(before): fail("rank_candidates changed existing columns")
if "region_z" not in o: fail("missing feature not reported")
yq = pd.read_csv(f"{B}/labels/consensus.csv"); yq = yq[yq.study_id == "Q06"].consensus_class.eq("cSS").values
top = after.sort_values("rank_model").head(int(yq.sum())).cand_id.values
if not yq[top - 1].all(): fail("model ranking on the test study does not put the cSS first")
if run(f"{S}/rank_candidates.py", "Q00", ok=False).returncode == 0: fail("a training study was ranked without --allow-training-study")
print("rank_candidates: columns appended at the end, others unchanged; held-out test study ranked cSS first; "
      "missing feature filled and named; training study refused")
# 2. no signal
make(False, 3, n_pat=8)
run(f"{S}/train_ranker.py", "--boot", "300")
M = json.load(open(f"{B}/review/ranker_model.json"))
lo, hi = M["lopo_auc_ci95"]
if not (lo <= 0.5 <= hi) or abs(M["lopo_auc"] - 0.5) > 0.15: fail(f"noise labels: LOPO AUC {M['lopo_auc']} {M['lopo_auc_ci95']}")
print(f"noise labels: LOPO AUC {M['lopo_auc']} {M['lopo_auc_ci95']} (interval includes 0.5)")
# too few patients
make(True, 4, n_pat=2, test_pat=0)
r = run(f"{S}/train_ranker.py", ok=False)
if r.returncode == 0 or "need >= 3" not in (r.stdout + r.stderr): fail("2 patients were accepted")
print("2 labelled patients: refused (needs >= 3)")
print("RESULT ranker OK")
