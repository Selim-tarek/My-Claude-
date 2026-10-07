#!/usr/bin/env python3
"""Interactive cSS review: shows each suspected candidate, you agree or disagree.
usage: review_css.py SUBJECT [--top N] [--redo]
  default: ALL candidates (v4). With --top N, candidates below N are recorded as NOT reviewed and
  the score is reported with that caveat (v3 silently scored them as "not cSS").

Keys (or click the buttons):
  y = cSS (agree)      v = vein      o = normal cortex      a = artifact
  i = near ICH         u = unsure    b = back one            [ / ] = move slice down/up
  q = finish and score
Decisions are saved after every key press, so you can quit and resume later.
At the end, accepted candidates are scored automatically (score_css.py)."""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import sys, os, subprocess, numpy as np, nibabel as nib, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import load_calls, save_calls
import matplotlib
if os.environ.get("CSS_REVIEW_TEST"): matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.patches import Rectangle

subj = sys.argv[1]
top = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else None
redo = "--redo" in sys.argv
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
csv = f"{base}/review/{subj}_candidates.csv"
calls_csv = f"{base}/review/{subj}_calls.csv"

swi = nib.as_closest_canonical(nib.load(f"{base}/data/{subj}_swi.nii"))      # RAS for display
cnd = nib.as_closest_canonical(nib.load(f"{base}/review/{subj}_candidates.nii.gz"))
I = swi.get_fdata().astype(np.float32); C = cnd.get_fdata().astype(int)
vox = swi.header.get_zooms()[:3]
br = I > 0
lo, hi = np.percentile(I[br], [1, 99]) if br.any() else (I.min(), I.max())
df = pd.read_csv(csv)
ids = df.cand_id.astype(int).tolist()[:top]
calls = {}
if not redo:
    # v4 fix: calls are matched by candidate fingerprint (centroid + volume), not by rank, so a
    # re-run of the detector can no longer attach an old call to a different lesion
    calls, warn = load_calls(base, subj, df)
    if warn: print("WARNING:", warn)
KEYS = {"y": "cSS", "v": "Vein", "o": "Normal", "a": "Artifact", "i": "Near ICH", "u": "Unsure"}
COL = {"cSS": "#2e7d32", "Vein": "#1565c0", "Normal": "#6d6d6d", "Artifact": "#ef6c00",
       "Near ICH": "#8e24aa", "Unsure": "#c9a400"}
HW = int(round(30 / vox[0]))          # 60 mm zoom window

class Reviewer:
    def __init__(self):
        self.pos = next((n for n, c in enumerate(ids) if c not in calls), 0)
        self.shift = 0
        self.fig = plt.figure(figsize=(14, 8.2))
        gs = self.fig.add_gridspec(2, 4, left=0.02, right=0.98, top=0.88, bottom=0.14, wspace=0.05, hspace=0.12)
        self.ov = self.fig.add_subplot(gs[:, 0])
        self.ax = [[self.fig.add_subplot(gs[r, c + 1]) for c in range(3)] for r in range(2)]
        self.btns = []
        labels = [("cSS (y)", "y"), ("Vein (v)", "v"), ("Normal (o)", "o"), ("Artifact (a)", "a"),
                  ("Near ICH (i)", "i"), ("Unsure (u)", "u"), ("◀ Back (b)", "b"),
                  ("Slice ↓ [", "["), ("Slice ↑ ]", "]"), ("Finish (q)", "q")]
        w = 0.094
        for n, (lab, k) in enumerate(labels):
            bax = self.fig.add_axes([0.02 + n * (w + 0.003), 0.03, w, 0.06])
            b = Button(bax, lab, color=COL.get(KEYS.get(k, ""), "#e0e0e0") if k in KEYS else "#e0e0e0",
                       hovercolor="#ffffff")
            if k in KEYS: b.label.set_color("white"); b.label.set_fontweight("bold")
            b.on_clicked(lambda e, k=k: self.key(k)); self.btns.append(b)
        self.fig.canvas.mpl_connect("key_press_event", lambda e: self.key(e.key))
        self.show()

    def show(self):
        if self.pos >= len(ids): return self.finish()
        cid = ids[self.pos]; r = df[df.cand_id == cid].iloc[0]
        m = C == cid
        if not m.any():
            self.pos += 1; return self.show()
        pts = np.argwhere(m); ci, cj = pts[:, 0].mean(), pts[:, 1].mean()
        ks = np.bincount(pts[:, 2]); kc = int(ks.argmax()) + self.shift
        kc = int(np.clip(kc, 1, I.shape[2] - 2))
        i0, i1 = int(max(ci - HW, 0)), int(min(ci + HW, I.shape[0]))
        j0, j1 = int(max(cj - HW, 0)), int(min(cj + HW, I.shape[1]))
        # overview
        self.ov.clear(); self.ov.imshow(I[:, :, kc].T, cmap="gray", origin="lower", vmin=lo, vmax=hi)
        self.ov.add_patch(Rectangle((i0, j0), i1 - i0, j1 - j0, fill=False, ec="yellow", lw=1.5))
        self.ov.set_title("whole slice (yellow = zoom)", fontsize=9); self.ov.axis("off")
        self.ov.text(2, I.shape[1] - 6, "L", color="yellow", fontsize=11); self.ov.text(I.shape[0] - 12, I.shape[1] - 6, "R", color="yellow", fontsize=11)
        # zoomed: top row raw, bottom row with outline; slices kc-1, kc, kc+1
        for c, k in enumerate((kc - 1, kc, kc + 1)):
            for row in (0, 1):
                a = self.ax[row][c]; a.clear()
                a.imshow(I[i0:i1, j0:j1, k].T, cmap="gray", origin="lower", vmin=lo, vmax=hi)
                if row == 1 and m[i0:i1, j0:j1, k].any():
                    a.contour(m[i0:i1, j0:j1, k].T.astype(float), levels=[0.5], colors="red", linewidths=1.2)
                a.set_xticks([]); a.set_yticks([])
                if row == 0: a.set_title(f"slice {k}" + ("  (centre)" if k == kc else ""), fontsize=9)
        self.ax[0][0].set_ylabel("raw SWI", fontsize=9); self.ax[1][0].set_ylabel("suspect outlined", fontsize=9)
        flags = [n for n, f in (("artifact zone", "artifact_zone"), ("midline", "midline_zone"), ("vein-like", "vein_like")) if int(r.get(f, 0) or 0)]
        done = sum(1 for c in ids if c in calls)
        prev = calls.get(cid)
        self.fig.suptitle(f"{subj}   candidate #{cid}  ({self.pos + 1}/{len(ids)}, {done} decided)   —   "
                          f"{r.hemi} {r.region}   {r.volume_mm3} mm³, {int(r.n_slices)} slices, darkness z {r.darkness_z}"
                          + (f"   [{', '.join(flags)}]" if flags else "")
                          + (f"\nyour previous call: {prev}" if prev else "\nIs the red outline cSS?"),
                          fontsize=11, color=COL.get(prev, "black"))
        self.fig.canvas.draw_idle()

    def key(self, k):
        if k in KEYS:
            calls[ids[self.pos]] = KEYS[k]; self.save(); self.pos += 1; self.shift = 0; self.show()
        elif k == "b":
            self.pos = max(self.pos - 1, 0); self.shift = 0; self.show()
        elif k == "[": self.shift -= 1; self.show()
        elif k == "]": self.shift += 1; self.show()
        elif k == "q": self.finish()

    def save(self):
        save_calls(base, subj, df, calls)

    def finish(self):
        self.save(); plt.close(self.fig)

R = Reviewer()
if os.environ.get("CSS_REVIEW_TEST"):
    R.fig.savefig(f"{base}/review/{subj}_review_preview.png", dpi=80)
    for n, c in enumerate(ids): calls[c] = ["cSS", "Vein", "Near ICH", "Normal"][n % 4]
    R.save()
else:
    plt.show()

# ---- write decisions and score
accepted = sorted(c for c, v in calls.items() if v == "cSS")
ich = sorted(c for c, v in calls.items() if v == "Near ICH")
unsure = sorted(c for c, v in calls.items() if v == "Unsure")
print(f"\n{subj}: reviewed {len(calls)} of top {len(ids)}  |  cSS {accepted}  |  near-ICH {ich}"
      + (f"  |  unsure {unsure} (counted as NOT cSS - re-check with: review_css.py {subj})" if unsure else ""))
not_rev = len(df) - len(calls)
if not_rev:
    print(f"WARNING: {not_rev} of {len(df)} candidates NOT reviewed - the score is a lower bound "
          f"(continue with: review_css.py {subj})")
cmd = [sys.executable, f"{base}/scripts/mark_css.py", subj, ",".join(map(str, accepted)) or "none",
       "--reviewed", ",".join(map(str, sorted(calls))) or "none"]
if ich: cmd += ["--ich", ",".join(map(str, ich))]
subprocess.run(cmd, check=False)
print(f"decisions saved: {calls_csv}")
