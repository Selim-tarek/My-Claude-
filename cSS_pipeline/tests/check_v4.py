"""Check v4 detector output on the phantom_v4 truth.   usage: check_v4.py NAME
Prints each candidate's truth class and features, sensitivity, and how well score_v3 / score_v4
and each definition feature separate true cSS from mimics. Exits 1 if a sanity check fails."""
import os, sys, numpy as np, nibabel as nib, pandas as pd
B = os.environ["CSS_BASE"]; name = sys.argv[1]
ld = lambda p: nib.load(p).get_fdata().astype(int)
truth, veins = ld(f"{B}/work/{name}_truth.nii.gz"), ld(f"{B}/work/{name}_veins.nii.gz")
cand = ld(f"{B}/review/{name}_candidates.nii.gz")
df = pd.read_csv(f"{B}/review/{name}_candidates.csv")
cls = []
for c in df.cand_id:
    m = cand == c; t = np.bincount(truth[m], minlength=7)[1:]; v = np.bincount(veins[m], minlength=7)[1:]
    if t.sum() >= v.sum() and t.sum() > 0.2 * m.sum(): cls.append(f"cSS{t.argmax() + 1}")
    elif v.sum() > 0: cls.append(f"vein{v.argmax() + 1}")
    else: cls.append("other")
df["truth"] = cls
F = ["darkness_z", "pial_dist_mm", "bank_frac", "tube_ratio", "surface_alignment", "tram_frac",
     "vein_tree_mm", "mirror_dark_frac", "score_v3", "score_v4"]
print(df[["cand_id", "truth"] + F].to_string(index=False))
# a lesion is found when kept candidates cover >=30 % of it (merged neighbours count for both)
found = [t for t in range(1, int(truth.max()) + 1) if (cand[truth == t] > 0).mean() >= 0.3]
n = int(truth.max()); print(f"\nsensitivity: {len(found)}/{n}  found lesions {found}  (coverage >=30 %)")
pos = df[df.truth.str.startswith("cSS")]; neg = df[~df.truth.str.startswith("cSS")]
def auc(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float); x, y = x[~np.isnan(x)], y[~np.isnan(y)]
    if not len(x) or not len(y): return float("nan")
    return float((x[:, None] > y[None, :]).mean() + 0.5 * (x[:, None] == y[None, :]).mean())
print("AUC cSS vs mimics (>0.5 = higher in cSS):  " +
      "  ".join(f"{f} {auc(pos[f], neg[f]):.2f}" for f in F))
vk = [v for v in range(1, int(veins.max()) + 1) if (cand[veins == v] > 0).mean() >= 0.3]
print(f"veins still shown as candidates: {len(vk)}/{int(veins.max())} {vk}")
ich = nib.load(f"{B}/work/{name}_ich_used.nii.gz").get_fdata() > 0
icht = nib.load(f"{B}/work/{name}_ich_truth.nii.gz").get_fdata() > 0
ich_dice = 2 * (ich & icht).sum() / max(ich.sum() + icht.sum(), 1)
print(f"automatic ICH mask Dice vs truth: {ich_dice:.2f}")
ok = len(found) >= n - 1 and auc(pos.score_v4, neg.score_v4) >= 0.8 and ich_dice > 0.5
print("RESULT v4 " + ("OK" if ok else "FAIL") +
      f" found={len(found)} total={n} auc_v3={auc(pos.score_v3, neg.score_v3):.2f} "
      f"auc_v4={auc(pos.score_v4, neg.score_v4):.2f} ich_dice={ich_dice:.2f}")
sys.exit(0 if ok else 1)
