"""v5.1.3 checks: review-priority tiers (scripts/css_tiers.py) - a reading order, never an exclusion.
usage: check_v513.py    (units need nothing; the phantom part needs PH3 / PH6 aligned in $CSS_BASE)
1. every branch on hand-made rows: QC insufficient -> 5; mixed object -> 2; microbleed-like / artifact zone -> 4;
   calcification / FLAIR context -> 2; convergent vein evidence -> 3; a single vein-ish feature alone -> still 1.
2. PH3 / PH6: every true cSS lesion tier 1, every vein (shown or excluded) tier 3 - on these phantoms, whose veins
   were built to differ; not real-world performance. The tier changes no candidate or exclusion.
3. review_css.py --by-priority (headless test mode) orders by tier; the files it writes are restored afterwards."""
import os, sys, subprocess, numpy as np, nibabel as nib, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
from css_tiers import review_priority as T
def fail(m): print("FAIL:", m); sys.exit(1)
css = dict(tube_ratio=0.25, pial_dist_mm=-1.1, csf_depth_mm=0.0, local_cont_max_mm=1.0, terminal_ends=2,
           end_type_max="stops", end_type_min="stops", vessel_wm_mm=0.0, elongation_3d=3.0, mimic_class="none")
cases = [("css", css, "ok", 1), ("qc", css, "insufficient_quality", 5), ("lowconf keeps tiers", css, "low_confidence_review", 1),
         ("mixed", dict(css, mixed_object=1), "ok", 2), ("mb", dict(css, mimic_class="mb_cluster"), "ok", 4),
         ("cmb", dict(css, cmb_like=1), "ok", 4), ("artifact", dict(css, artifact_zone=1), "ok", 4),
         ("calcium", dict(css, mimic_class="calcification"), "ok", 2), ("flair", dict(css, mimic_class="acute"), "ok", 2),
         ("transcortical", dict(css, vessel_wm_mm=6.0, local_cont_max_mm=9.0), "ok", 3),
         ("sinus", dict(css, end_type_max="midline", local_cont_max_mm=12.0), "ok", 3),
         ("through-running tube", dict(css, tube_ratio=0.5, terminal_ends=0), "ok", 3),
         ("off-cortex vessel", dict(css, pial_dist_mm=1.2, csf_depth_mm=1.0, elongation_3d=20.0), "ok", 3),
         ("one feature: off the cortex only", dict(css, pial_dist_mm=1.2), "ok", 1),
         ("one feature: tubular only", dict(css, tube_ratio=0.6), "ok", 1),
         ("one feature: WM continuation, short and flat", dict(css, vessel_wm_mm=6.0), "ok", 1)]
for name, row, qc, want in cases:
    got = T(row, qc)[0]
    if got != want: fail(f"{name}: tier {got}, expected {want} ({T(row, qc)[2]})")
print(f"units: {len(cases)} hand-made rows - each branch as specified; a single vein-ish feature alone stays tier 1")

B = os.environ.get("CSS_BASE")
if not B or not os.path.exists(f"{B}/data/PH3_swi.nii"):
    print("phantom part skipped"); print("RESULT v513 OK"); sys.exit(0)
S = f"{B}/scripts"; py = sys.executable
ve = nib.load(f"{B}/work/PH3_veins.nii.gz").get_fdata().astype(int); seen = []
for s in ("PH3", "PH6"):
    if not os.path.exists(f"{B}/work/{s}_seg_swispace.nii.gz"): continue
    r = subprocess.run([py, f"{S}/detect_css.py", s], capture_output=True, text=True)
    if r.returncode: fail(f"detect_css {s}")
    tr = nib.load(f"{B}/work/{s}_truth.nii.gz").get_fdata().astype(int)
    sv = nib.load(f"{B}/work/PH6_surfvein.nii.gz").get_fdata().astype(int) if s == "PH6" else np.zeros_like(tr)
    for lst, idc in (("candidates", "cand_id"), ("excluded", "excl_id")):
        d = pd.read_csv(f"{B}/review/{s}_{lst}.csv"); M = nib.load(f"{B}/review/{s}_{lst}.nii.gz").get_fdata().astype(int)
        for _, rr in d.iterrows():
            m = M == rr[idc]
            if (tr[m] > 0).mean() > 0.2 and rr.review_priority_rank != 1: fail(f"{s}: a true cSS lesion in tier {rr.review_priority_rank}")
            if ((ve[m] > 0).mean() > 0.3 or (sv[m] > 0).mean() > 0.5) and (tr[m] > 0).mean() <= 0.2 and rr.review_priority_rank != 3:
                fail(f"{s}: a phantom vein in tier {rr.review_priority_rank} ({rr.review_priority_reason})")
    seen.append(s)
print(f"phantoms {', '.join(seen)}: every cSS lesion tier 1, every vein tier 3 (constructed veins; not performance)")
env = dict(os.environ, CSS_REVIEW_TEST="1")
# the headless reviewer writes calls, a mask and a score: keep the shared test folder as it was for later tests
touched = [f"{B}/review/{f}" for f in ("PH3_calls.csv", "PH3_candidates.csv", "PH3_score.json", "PH3_css_mask.nii.gz",
                                       "PH3_review_preview.png", "css_scores.csv")]
saved = {f: open(f, "rb").read() for f in touched if os.path.exists(f)}
r = subprocess.run([py, f"{S}/review_css.py", "PH3", "--by-priority", "--redo", "--top", "3"], capture_output=True, text=True, env=env)
for f in touched:
    if f in saved: open(f, "wb").write(saved[f])
    elif os.path.exists(f): os.remove(f)
if "order: review priority" not in r.stdout: fail(f"--by-priority not applied:\n{r.stdout[-500:]}{r.stderr[-500:]}")
print("review_css.py --by-priority: order by tier applied (" + [l for l in r.stdout.splitlines() if l.startswith("order:")][0] + ")")
print("RESULT v513 OK")
