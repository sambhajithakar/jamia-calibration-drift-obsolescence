"""
Fit frozen models in the base period, follow them forward across evaluation cycles
under five maintenance policies, and record design-based metrics per cycle.

Executes PRESPECIFICATION_v1.md (sha256 5daf3197b26be2714c3f5f2b6fc240873792ce5f487360869f677e598af4aac8).
"""
import os, json, warnings
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, log_loss

warnings.filterwarnings("ignore")
RNG = np.random.default_rng(20261106)
OUT = "/home/claude/jamia2026/03_analysis/out"

CYCLES = ["1999-2000","2001-2002","2003-2004","2005-2006","2007-2008","2009-2010",
          "2011-2012","2013-2014","2015-2016","2017-2018","2021-2023"]
MIDYEAR = {c: y for c, y in zip(CYCLES, [1999.5,2001.5,2003.5,2005.5,2007.5,2009.5,
                                         2011.5,2013.5,2015.5,2017.5,2022.0])}
BASE = CYCLES[:2]
EVAL = CYCLES[2:]

Q   = ["age","female","race","educ_hs","pir","pir_miss","smoke_ever","smoke_now"]
QE  = Q + ["bmi","waist","height","weight"]
TASKS = ["y_undiag_diabetes","y_obesity","y_central_obesity","y_htn_measured",
         "y_arthritis","y_chf","y_chd","y_mi","y_stroke","y_cancer"]
# tasks whose outcome is derived from the examination block: Q only, to avoid using the outcome
Q_ONLY_TASKS = {"y_obesity","y_central_obesity"}
# blood pressure method changes in 2021-2023 (pre-specified)
STOP_2017 = {"y_htn_measured"}

df = pd.read_csv(f"{OUT}/analytic.csv.gz")
df["pir_miss"] = df.pir.isna().astype(float)

# impute once within cycle: median continuous, mode binary (pre-specified)
CONT = ["age","pir","bmi","waist","height","weight"]
BIN  = ["female","educ_hs","smoke_ever","smoke_now","race"]
for c in CONT:
    df[c] = df.groupby("cycle")[c].transform(lambda s: s.fillna(s.median()))
for c in BIN:
    df[c] = df.groupby("cycle")[c].transform(
        lambda s: s.fillna(s.mode().iloc[0] if not s.mode().empty else 0))

def wmean(x, w):
    return float(np.sum(w * x) / np.sum(w))

def wauc(y, p, w):
    # weighted AUROC via the Mann-Whitney form
    o = np.argsort(p); y, p, w = y[o], p[o], w[o]
    pos, neg = y == 1, y == 0
    if pos.sum() == 0 or neg.sum() == 0:
        return np.nan
    try:
        return float(roc_auc_score(y, p, sample_weight=w))
    except Exception:
        return np.nan

def net_benefit(y, p, w, t):
    """standardised-scale net benefit at threshold t, and for treat-all"""
    W = np.sum(w)
    flag = p >= t
    tp = np.sum(w[(flag) & (y == 1)]) / W
    fp = np.sum(w[(flag) & (y == 0)]) / W
    nb = tp - fp * t / (1 - t)
    prev = np.sum(w[y == 1]) / W
    nb_all = prev - (1 - prev) * t / (1 - t)
    return nb, nb_all, prev

def calib(y, p, w):
    """calibration-in-the-large and slope from a weighted logistic recalibration"""
    eps = 1e-6
    lp = np.log(np.clip(p, eps, 1 - eps) / (1 - np.clip(p, eps, 1 - eps)))
    try:
        m = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000)
        m.fit(lp.reshape(-1, 1), y, sample_weight=w)
        slope = float(m.coef_[0][0]); inter = float(m.intercept_[0])
    except Exception:
        return np.nan, np.nan, None
    # calibration-in-the-large: mean predicted vs observed on the log-odds scale
    citl_m = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000)
    try:
        citl_m.fit(np.zeros((len(y), 1)), y, sample_weight=w)
    except Exception:
        pass
    obs = np.sum(w * y) / np.sum(w)
    exp = np.sum(w * p) / np.sum(w)
    citl = float(np.log(max(obs, eps) / max(1 - obs, eps)) -
                 np.log(max(exp, eps) / max(1 - exp, eps)))
    return citl, slope, (inter, slope)

def calib_intercept_only(y, p, w):
    """intercept-only recalibration: fit a in  y ~ sigmoid(a + logit(p)), slope fixed at 1"""
    from scipy.optimize import minimize_scalar
    eps = 1e-6
    lp = np.log(np.clip(p, eps, 1-eps) / (1 - np.clip(p, eps, 1-eps)))
    def nll(a):
        z = a + lp
        return float(np.sum(w * (np.logaddexp(0, z) - y * z)))
    r = minimize_scalar(nll, bounds=(-10, 10), method="bounded")
    return float(r.x)


def apply_recal(p, inter, slope):
    eps = 1e-6
    lp = np.log(np.clip(p, eps, 1 - eps) / (1 - np.clip(p, eps, 1 - eps)))
    z = inter + slope * lp
    return 1 / (1 + np.exp(-z))

def fit(X, y, w, family):
    if family == "lr":
        m = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=2000)
    else:
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.06,
                                           max_leaf_nodes=15, random_state=0)
    m.fit(X, y, sample_weight=w)
    return m

rows, elig_report = [], []
for task in TASKS:
    base = df[df.cycle.isin(BASE) & df[task].notna()]
    n_ev = int(base[task].sum()); prev = float(base[task].mean())
    eligible = (prev >= 0.01) and (n_ev >= 100)
    elig_report.append(dict(task=task, base_n=len(base), base_events=n_ev,
                            base_prev=round(prev, 4), eligible=eligible))
    if not eligible:
        continue

    thr = max(round(prev, 2), 0.01)        # pre-specified: base-period prevalence
    predsets = {"Q": Q} if task in Q_ONLY_TASKS else {"Q": Q, "QE": QE}
    evals = [c for c in EVAL if not (task in STOP_2017 and c == "2021-2023")]

    for ps_name, feats in predsets.items():
        Xb, yb = base[feats].to_numpy(float), base[task].to_numpy(float)
        wb = base.WTMEC2YR.to_numpy(float); wb = wb / wb.mean()

        for family in ["lr", "gbt"]:
            model = fit(Xb, yb, wb, family)
            state = {p: None for p in ["P1", "P2", "P3", "P4"]}   # recal params / refit model
            hist = []

            for cyc in evals:
                ev = df[(df.cycle == cyc) & df[task].notna()]
                if len(ev) < 200:
                    continue
                Xe, ye = ev[feats].to_numpy(float), ev[task].to_numpy(float)
                we = ev.WTMEC2YR.to_numpy(float); we = we / we.mean()
                psu, stra = ev.SDMVPSU.to_numpy(), ev.SDMVSTRA.to_numpy()

                p0 = model.predict_proba(Xe)[:, 1]
                preds = {"P0": p0}
                preds["P1"] = apply_recal(p0, state["P1"][0], 1.0) if state["P1"] else p0
                preds["P2"] = apply_recal(p0, *state["P2"]) if state["P2"] else p0
                preds["P3"] = apply_recal(p0, *state["P3"]) if state["P3"] else p0
                preds["P4"] = (state["P4"].predict_proba(Xe)[:, 1]
                               if state["P4"] is not None else p0)

                for pol, pp in preds.items():
                    auc = wauc(ye, pp, we)
                    citl, slope, _ = calib(ye, pp, we)
                    nb, nb_all, prv = net_benefit(ye, pp, we, thr)
                    rows.append(dict(task=task, predset=ps_name, family=family,
                                     policy=pol, cycle=cyc, midyear=MIDYEAR[cyc],
                                     n=len(ev), prev=prv, thr=thr, auroc=auc,
                                     citl=citl, slope=slope, nb=nb, nb_all=nb_all,
                                     nb_default=max(nb_all, 0.0),
                                     fail=bool(nb <= max(nb_all, 0.0))))

                # update maintenance state using THIS cycle, applied from the next one
                citl0, slope0, par0 = calib(ye, p0, we)
                if par0 is not None:
                    state["P1"] = (calib_intercept_only(ye, p0, we), 1.0)
                    state["P2"] = par0
                    if not (-0.2 <= citl0 <= 0.2):
                        state["P3"] = par0
                state["P4"] = fit(Xe, ye, we, family)
                hist.append(cyc)

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/followup_metrics.csv", index=False)
pd.DataFrame(elig_report).to_csv(f"{OUT}/task_eligibility.csv", index=False)

print("ELIGIBILITY")
print(pd.DataFrame(elig_report).to_string(index=False))
print()
print(f"model-task combinations: {res.groupby(['task','predset','family']).ngroups}")
print(f"rows of follow-up metrics: {len(res):,}")
