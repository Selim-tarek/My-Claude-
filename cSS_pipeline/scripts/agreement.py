#!/usr/bin/env python3
"""Reproducibility statistics as in van Harten et al. (NeuroImage Clin 2023;38:103447).
usage: agreement.py --a DIR_A --b DIR_B [--col grown_volume_mm3] [--plot out.png]
       agreement.py --vs-score DIR [--col surface_mm2] [--plot out.png]
  DIR = a copy of the review folder after one complete scoring session, e.g.
        cp -R ~/css_project/review ~/css_project/review_reader1_session1
        (it must contain css_scores.csv and, for Dice, ID_css_mask.nii.gz)
  --a/--b      two sessions of one reader (intra-observer) or two readers (inter-observer):
               Pearson r, ICC(A,1) absolute agreement (two-way), Bland-Altman bias and 95 % limits of
               agreement, and per-subject Dice of the cSS masks. van Harten: ICC 0.995, Pearson 0.991,
               mean Dice 0.75 (Dice is low for small, sparse masks even when volumes agree).
  --vs-score   size measure against the 0-4 multifocality score (their Fig. 4: wide spread inside the
               top categories = the ceiling effect a continuous measure avoids)
Compare sizes only between scans acquired with the same protocol (field strength, TE, voxel size)."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, numpy as np, pandas as pd, nibabel as nib

arg = lambda f, d="": sys.argv[sys.argv.index(f) + 1] if f in sys.argv else d
col = arg("--col", "grown_volume_mm3")
ex = lambda p: os.path.expanduser(p)


def icc_a1(x, y):
    """two-way, absolute agreement, single measures (McGraw & Wong ICC(A,1))"""
    Y = np.c_[x, y]; n, k = Y.shape
    gm = Y.mean(); rm = Y.mean(1); cm = Y.mean(0)
    msr = k * ((rm - gm) ** 2).sum() / (n - 1)
    msc = n * ((cm - gm) ** 2).sum() / (k - 1)
    mse = ((Y - rm[:, None] - cm[None, :] + gm) ** 2).sum() / ((n - 1) * (k - 1))
    return float((msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n))


if arg("--vs-score"):
    d = pd.read_csv(f"{ex(arg('--vs-score'))}/css_scores.csv")
    if col not in d.columns: sys.exit(f"column {col} not in css_scores.csv (re-score with score_css.py v4.13)")
    print(f"{col} by multifocality score ({len(d)} subjects):")
    for s, g in d.groupby("multifocality_0_4"):
        print(f"  score {s}: n={len(g)}  median {g[col].median():.1f}  range {g[col].min():.1f}-{g[col].max():.1f}")
    if arg("--plot"):
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 4))
        jit = np.random.default_rng(0).uniform(-0.12, 0.12, len(d))
        ax.scatter(d.multifocality_0_4 + jit, d[col], s=25)
        ax.set_xlabel("cSS multifocality score (0-4)"); ax.set_ylabel(col); ax.set_xticks(range(5))
        fig.tight_layout(); fig.savefig(ex(arg("--plot")), dpi=150); print(f"plot -> {arg('--plot')}")
    sys.exit(0)

if not (arg("--a") and arg("--b")): sys.exit(__doc__)
A, Bd = ex(arg("--a")), ex(arg("--b"))
da, db = pd.read_csv(f"{A}/css_scores.csv"), pd.read_csv(f"{Bd}/css_scores.csv")
m = da[["subject", col]].merge(db[["subject", col]], on="subject", suffixes=("_a", "_b")).dropna()
if len(m) < 2: sys.exit(f"need >=2 subjects scored in both folders (found {len(m)})")
x, y = m[f"{col}_a"].to_numpy(float), m[f"{col}_b"].to_numpy(float)
diff = x - y; bias = diff.mean(); sd = diff.std(ddof=1)
r = float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else float("nan")
print(f"{len(m)} subjects, {col}:  Pearson r {r:.3f}   ICC(A,1) {icc_a1(x, y):.3f}   "
      f"Bland-Altman bias {bias:.1f}, 95 % limits {bias - 1.96 * sd:.1f} to {bias + 1.96 * sd:.1f}")
dice = []
for s in m.subject:
    pa, pb = f"{A}/{s}_css_mask.nii.gz", f"{Bd}/{s}_css_mask.nii.gz"
    if os.path.exists(pa) and os.path.exists(pb):
        a_, b_ = nib.load(pa).get_fdata() > 0, nib.load(pb).get_fdata() > 0
        if a_.shape == b_.shape and (a_.any() or b_.any()):
            dice.append((s, 2 * (a_ & b_).sum() / (a_.sum() + b_.sum())))
if dice:
    dv = np.array([v for _, v in dice])
    print(f"Dice of the cSS masks: mean {dv.mean():.2f} +/- {dv.std():.2f}  ("
          + ", ".join(f"{s} {v:.2f}" for s, v in dice) + ")")
print(m.to_string(index=False))
if arg("--plot"):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    lim = [0, max(x.max(), y.max()) * 1.05]
    ax[0].scatter(x, y); ax[0].plot(lim, lim, "k--", lw=0.8)
    ax[0].set_xlabel(f"{col} (A)"); ax[0].set_ylabel(f"{col} (B)"); ax[0].set_title(f"r = {r:.3f}")
    mean = (x + y) / 2; ax[1].scatter(mean, diff)
    for v, ls in ((bias, "-"), (bias - 1.96 * sd, "--"), (bias + 1.96 * sd, "--")): ax[1].axhline(v, color="k", ls=ls, lw=0.8)
    ax[1].set_xlabel("mean of A and B"); ax[1].set_ylabel("A - B"); ax[1].set_title("Bland-Altman")
    fig.tight_layout(); fig.savefig(ex(arg("--plot")), dpi=150); print(f"plot -> {arg('--plot')}")
