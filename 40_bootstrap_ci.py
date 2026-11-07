"""
Design-based confidence intervals for the frozen-model (P0) metrics.

Executes PRESPECIFICATION_v1.md section 8: "Confidence intervals by bootstrap resampling
of PSUs within strata, 500 replicates, seed fixed at 20261106."

The per-cell seed is derived from SEED and a SHA-256 digest of the cell key, so the
replicate draws are identical on every machine and in every process.

Models are refit from the base period exactly as in 10_fit_and_follow.py (same code, same
seeds), so the point estimates written here must reproduce out/followup_metrics.csv for
policy P0. That equality is asserted at the end of the run.

Writes out/bootstrap_ci.csv — one row per task/predset/family/cycle with the point
estimate and the 2.5th/97.5th percentile of the replicate distribution for AUROC,
calibration-in-the-large, calibration slope and net benefit.
"""
import hashlib, warnings
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
OUT = "/home/claude/jamia2026/03_analysis/out"
NBOOT = 500
SEED = 20261106

CYCLES = ["1999-2000","2001-2002","2003-2004","2005-2006","2007-2008","2009-2010",
          "2011-2012","2013-2014","2015-2016","2017-2018","2021-2023"]
MIDYEAR = {c: y for c, y in zip(CYCLES, [1999.5,2001.5,2003.5,2005.5,2007.5,2009.5,
                                         2011.5,2013.5,2015.5,2017.5,2022.0])}
BASE, EVAL = CYCLES[:2], CYCLES[2:]

Q  = ["age","female","race","educ_hs","pir","pir_miss","smoke_ever","smoke_now"]
QE = Q + ["bmi","waist","height","weight"]
TASKS = ["y_undiag_diabetes","y_obesity","y_central_obesity","y_htn_measured",
         "y_arthritis","y_chf","y_chd","y_mi","y_stroke","y_cancer"]
Q_ONLY_TASKS = {"y_obesity","y_central_obesity"}
STOP_2017 = {"y_htn_measured"}

df = pd.read_csv(f"{OUT}/analytic.csv.gz")
df["pir_miss"] = df.pir.isna().astype(float)
CONT = ["age","pir","bmi","waist","height","weight"]
BIN  = ["female","educ_hs","smoke_ever","smoke_now","race"]
for c in CONT:
    df[c] = df.groupby("cycle")[c].transform(lambda s: s.fillna(s.median()))
for c in BIN:
    df[c] = df.groupby("cycle")[c].transform(
        lambda s: s.fillna(s.mode().iloc[0] if not s.mode().empty else 0))


def wauc(y, p, w):
    if (y == 1).sum() == 0 or (y == 0).sum() == 0:
        return np.nan
    try:
        return float(roc_auc_score(y, p, sample_weight=w))
    except Exception:
        return np.nan


def net_benefit(y, p, w, t):
    W = np.sum(w)
    flag = p >= t
    tp = np.sum(w[flag & (y == 1)]) / W
    fp = np.sum(w[flag & (y == 0)]) / W
    return tp - fp * t / (1 - t)


def citl_slope(y, p, w):
    eps = 1e-6
    pc = np.clip(p, eps, 1 - eps)
    lp = np.log(pc / (1 - pc))
    try:
        m = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000)
        m.fit(lp.reshape(-1, 1), y, sample_weight=w)
        slope = float(m.coef_[0][0])
    except Exception:
        slope = np.nan
    obs = np.sum(w * y) / np.sum(w)
    exp = np.sum(w * p) / np.sum(w)
    citl = float(np.log(max(obs, eps) / max(1 - obs, eps)) -
                 np.log(max(exp, eps) / max(1 - exp, eps)))
    return citl, slope


def fit(X, y, w, family):
    if family == "lr":
        m = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=2000)
    else:
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.06,
                                           max_leaf_nodes=15, random_state=0)
    m.fit(X, y, sample_weight=w)
    return m


def psu_bootstrap_indices(stra, psu, rng, nboot):
    """Resample PSUs with replacement within each stratum; return a list of row-index arrays.

    This is the standard with-replacement PSU bootstrap for a stratified multistage design.
    Strata contributing a single PSU cannot be resampled and are carried through unchanged,
    which is the conventional (slightly conservative) treatment.
    """
    idx = np.arange(len(stra))
    groups = {}
    for s in np.unique(stra):
        sel = idx[stra == s]
        units = np.unique(psu[sel])
        groups[s] = (units, {u: sel[psu[sel] == u] for u in units})
    out = []
    for _ in range(nboot):
        parts = []
        for s, (units, members) in groups.items():
            if len(units) == 1:
                parts.append(members[units[0]])
                continue
            drawn = rng.choice(units, size=len(units), replace=True)
            for u in drawn:
                parts.append(members[u])
        out.append(np.concatenate(parts))
    return out


rows = []
for task in TASKS:
    base = df[df.cycle.isin(BASE) & df[task].notna()]
    n_ev = int(base[task].sum()); prev = float(base[task].mean())
    if not ((prev >= 0.01) and (n_ev >= 100)):
        continue
    thr = max(round(prev, 2), 0.01)
    predsets = {"Q": Q} if task in Q_ONLY_TASKS else {"Q": Q, "QE": QE}
    evals = [c for c in EVAL if not (task in STOP_2017 and c == "2021-2023")]

    for ps_name, feats in predsets.items():
        Xb, yb = base[feats].to_numpy(float), base[task].to_numpy(float)
        wb = base.WTMEC2YR.to_numpy(float); wb = wb / wb.mean()
        for family in ["lr", "gbt"]:
            model = fit(Xb, yb, wb, family)
            for cyc in evals:
                ev = df[(df.cycle == cyc) & df[task].notna()]
                if len(ev) < 200:
                    continue
                Xe = ev[feats].to_numpy(float)
                ye = ev[task].to_numpy(float)
                we = ev.WTMEC2YR.to_numpy(float); we = we / we.mean()
                stra = ev.SDMVSTRA.to_numpy()
                psu = ev.SDMVPSU.to_numpy()
                p0 = model.predict_proba(Xe)[:, 1]

                pt_auc = wauc(ye, p0, we)
                pt_citl, pt_slope = citl_slope(ye, p0, we)
                pt_nb = net_benefit(ye, p0, we, thr)

                # A stable, process-independent per-cell seed. Python's hash() of a
                # tuple of strings is salted per process unless PYTHONHASHSEED is set,
                # so the previous version was not in fact reproducible across runs.
                key = "|".join((task, ps_name, family, cyc)).encode("utf8")
                rng = np.random.default_rng(
                    SEED + int(hashlib.sha256(key).hexdigest()[:12], 16) % 10_000_019)
                reps = psu_bootstrap_indices(stra, psu, rng, NBOOT)
                b_auc, b_citl, b_slope, b_nb = [], [], [], []
                for ix in reps:
                    yy, pp, ww = ye[ix], p0[ix], we[ix]
                    b_auc.append(wauc(yy, pp, ww))
                    c, s = citl_slope(yy, pp, ww)
                    b_citl.append(c); b_slope.append(s)
                    b_nb.append(net_benefit(yy, pp, ww, thr))

                def q(a):
                    a = np.asarray(a, float); a = a[np.isfinite(a)]
                    if a.size < 50:
                        return (np.nan, np.nan)
                    return (float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5)))

                la, ua = q(b_auc); lc, uc = q(b_citl)
                ls, us = q(b_slope); ln, un = q(b_nb)
                rows.append(dict(task=task, predset=ps_name, family=family,
                                 cycle=cyc, midyear=MIDYEAR[cyc],
                                 auroc=pt_auc, auroc_lo=la, auroc_hi=ua,
                                 citl=pt_citl, citl_lo=lc, citl_hi=uc,
                                 slope=pt_slope, slope_lo=ls, slope_hi=us,
                                 nb=pt_nb, nb_lo=ln, nb_hi=un))
            print(f"  done {task} {ps_name} {family}", flush=True)

out = pd.DataFrame(rows)
out.to_csv(f"{OUT}/bootstrap_ci.csv", index=False)
print(f"\nwrote bootstrap_ci.csv — {len(out)} evaluation points, {NBOOT} replicates each")

# ---- assert the point estimates reproduce the primary run exactly ----
prim = pd.read_csv(f"{OUT}/followup_metrics.csv")
prim = prim[prim.policy == "P0"]
m = out.merge(prim, on=["task", "predset", "family", "cycle"], suffixes=("", "_p"))
assert len(m) == len(out), f"merge lost rows: {len(m)} vs {len(out)}"
for c in ["auroc", "citl", "slope", "nb"]:
    d = np.nanmax(np.abs(m[c] - m[c + "_p"]))
    print(f"  max |{c} - primary| = {d:.3e}")
    assert d < 1e-9, f"{c} does not reproduce the primary run"
print("point estimates reproduce out/followup_metrics.csv exactly")

# ---- how many breaches are distinguishable from sampling error ----
out["slope_breach"] = ~out.slope.between(0.80, 1.25)
out["slope_breach_ci"] = (out.slope_hi < 0.80) | (out.slope_lo > 1.25)
out["citl_breach"] = ~out.citl.between(-0.20, 0.20)
out["citl_breach_ci"] = (out.citl_hi < -0.20) | (out.citl_lo > 0.20)
print()
print(f"evaluation points                     : {len(out)}")
print(f"slope outside [0.80, 1.25]            : {out.slope_breach.sum()} "
      f"({100*out.slope_breach.mean():.1f}%)")
print(f"  ... with the whole 95% CI outside   : {out.slope_breach_ci.sum()} "
      f"({100*out.slope_breach_ci.mean():.1f}%)")
print(f"CITL outside [-0.20, 0.20]            : {out.citl_breach.sum()} "
      f"({100*out.citl_breach.mean():.1f}%)")
print(f"  ... with the whole 95% CI outside   : {out.citl_breach_ci.sum()} "
      f"({100*out.citl_breach_ci.mean():.1f}%)")
