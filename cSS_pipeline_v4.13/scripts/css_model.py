"""Shared model code for the learned rankers (v5.2): class-balanced L2 logistic regression, patient-grouped
leave-one-patient-out (LOPO) with optional nested selection of the L2 strength, metrics, model hashing.
numpy / scipy only. Used by train_ranker.py and classify_surface_patches.py; tested in tests/check_v52.py.
Every split is by PATIENT (study_id): all objects of one study are always on the same side."""
import hashlib, json, numpy as np
from scipy.optimize import minimize

L2_GRID = (0.1, 1.0, 10.0)


def prep(X, med=None, mu=None, sd=None):
    """median imputation + standardisation; parameters from X (the training fold) unless given"""
    X = np.asarray(X, float)
    if med is None: med = np.nanmedian(X, axis=0); med = np.where(np.isfinite(med), med, 0.0)
    X = np.where(np.isfinite(X), X, med)
    if mu is None: mu = X.mean(0); sd = X.std(0); sd = np.where(sd > 1e-9, sd, 1.0)
    return (X - mu) / sd, med, mu, sd


def fit_lr(Z, y, l2):
    """class-balanced L2 logistic regression on standardised Z -> (w, b); raises if only one class"""
    y = np.asarray(y, float)
    if len(np.unique(y)) < 2: raise ValueError("both classes are required")
    cw = np.where(y == 1, 0.5 / y.mean(), 0.5 / (1 - y.mean()))
    def f(p):
        w, b = p[:-1], p[-1]; t = np.clip(Z @ w + b, -40, 40)
        g = cw * (1 / (1 + np.exp(-t)) - y)
        return (cw * (np.logaddexp(0, t) - y * t)).sum() / len(y) + 0.5 * l2 * (w @ w) / len(y), \
            np.r_[Z.T @ g / len(y) + l2 * w / len(y), g.sum() / len(y)]
    r = minimize(f, np.zeros(Z.shape[1] + 1), jac=True, method="L-BFGS-B")
    return r.x[:-1], float(r.x[-1])


def predict(X, med, mu, sd, w, b):
    Z = prep(X, med, mu, sd)[0]
    return 1 / (1 + np.exp(-(Z @ w + b)))


def auc(y, s):
    y = np.asarray(y); s = np.asarray(s, float); pos, neg = s[y == 1], s[y == 0]
    if not len(pos) or not len(neg): return np.nan
    return float(((pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()) / (len(pos) * len(neg)))


def average_precision(y, s):
    """PR-AUC as average precision (precision at each positive, ranked by score)"""
    y = np.asarray(y); s = np.asarray(s, float)
    if not (y == 1).any(): return np.nan
    o = np.argsort(-s, kind="stable"); yy = y[o]; tp = np.cumsum(yy)
    return float((tp[yy == 1] / (np.nonzero(yy == 1)[0] + 1)).mean())


def brier(y, p):
    return float(np.mean((np.asarray(p, float) - np.asarray(y, float)) ** 2))


def lopo(X, y, groups, l2=1.0, nested=False):
    """leave-one-patient-out predictions. nested=True picks l2 from L2_GRID inside each training fold by an inner LOPO
    (AUC). Folds whose training part has one class only are SKIPPED and reported (never fabricated).
    -> (out-of-fold scores, fold log)"""
    X = np.asarray(X, float); y = np.asarray(y); g = np.asarray(groups)
    out = np.full(len(y), np.nan); log = []
    for p in sorted(set(g)):
        te, tr = g == p, g != p
        if len(np.unique(y[tr])) < 2:
            log.append(dict(study_id=str(p), status="skipped", reason="training fold has one class")); continue
        lam = l2
        if nested:
            best = (-1, l2)
            for cand in L2_GRID:
                s_in, _ = lopo(X[tr], y[tr], g[tr], cand, nested=False)
                a = auc(y[tr][np.isfinite(s_in)], s_in[np.isfinite(s_in)])
                if np.isfinite(a) and a > best[0]: best = (a, cand)
            lam = best[1]
        Z, med, mu, sd = prep(X[tr]); w, b = fit_lr(Z, y[tr], lam)
        out[te] = predict(X[te], med, mu, sd, w, b)
        log.append(dict(study_id=str(p), status="evaluated", n_test=int(te.sum()), l2=lam))
    return out, log


def patient_bootstrap(df, fn, n=500, seed=1):
    """95 % interval of fn(sub-dataframe) over bootstrap resamples of PATIENTS (column study_id)"""
    rng = np.random.default_rng(seed); pats = df.study_id.unique(); vals = []
    import pandas as pd
    for _ in range(n):
        d = pd.concat([df[df.study_id == p] for p in rng.choice(pats, len(pats), replace=True)])
        v = fn(d)
        if np.isfinite(v): vals.append(v)
    return [round(float(np.percentile(vals, 2.5)), 3), round(float(np.percentile(vals, 97.5)), 3)] if vals else [None, None]


def sens_at_fp(df, score, fp_per_scan=1.0):
    """lesion sensitivity when, in every study, candidates are read down the score order until fp_per_scan false
    positives have been passed (the reading-burden view; needs columns study_id, y and the score)"""
    found = tot = 0
    for _, d in df.groupby("study_id"):
        d = d.sort_values(score, ascending=False); fp = np.cumsum(d.y.values == 0)
        keep = fp <= fp_per_scan; found += int(d.y.values[keep].sum()); tot += int(d.y.sum())
    return found / tot if tot else np.nan


def model_hash(model):
    """sha256 of the model content (weights, features, preprocessing), excluding the hash field itself"""
    m = {k: v for k, v in model.items() if k not in ("model_hash",)}
    return hashlib.sha256(json.dumps(m, sort_keys=True, default=str).encode()).hexdigest()
