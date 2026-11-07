"""
Three analyses that complete the pre-specification and test the paper's mechanistic claim.

A. PRESPECIFICATION_v1 section 4, not previously reported: sensitivity at 0.5x and 2x the
   primary decision threshold.
B. PRESPECIFICATION_v1 section 9, not previously reported: measured hypertension with the
   2021-2023 cycle added, despite the change from auscultatory to oscillometric
   measurement, and the difference against the pre-specified stop at 2017-2018.
C. A decomposition that tests the manuscript's mechanistic claim directly. The claim is
   that the hypertension failure came from a threshold left fixed while prevalence moved,
   not from the fitted function. If so, holding the frozen model's predictions exactly as
   they are and letting only the threshold track contemporaneous design-based prevalence
   should remove the loss. If the loss survives, the claim is wrong.

Writes out/threshold_mechanism.csv and prints a summary. Nothing here is used to change
any primary result; it is reported as confirmatory and clearly labelled post-hoc, except
A and B, which execute the pre-specification.
"""
import warnings
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier

warnings.filterwarnings("ignore")
OUT = "/home/claude/jamia2026/03_analysis/out"

CYCLES = ["1999-2000","2001-2002","2003-2004","2005-2006","2007-2008","2009-2010",
          "2011-2012","2013-2014","2015-2016","2017-2018","2021-2023"]
MIDYEAR = {c: y for c, y in zip(CYCLES, [1999.5,2001.5,2003.5,2005.5,2007.5,2009.5,
                                         2011.5,2013.5,2015.5,2017.5,2022.0])}
BASE, EVAL = CYCLES[:2], CYCLES[2:]
Q  = ["age","female","race","educ_hs","pir","pir_miss","smoke_ever","smoke_now"]
QE = Q + ["bmi","waist","height","weight"]
TASKS = ["y_undiag_diabetes","y_obesity","y_central_obesity","y_htn_measured",
         "y_arthritis","y_chf","y_chd","y_mi","y_stroke","y_cancer"]
Q_ONLY = {"y_obesity","y_central_obesity"}
STOP_2017 = {"y_htn_measured"}

df = pd.read_csv(f"{OUT}/analytic.csv.gz")
df["pir_miss"] = df.pir.isna().astype(float)
for c in ["age","pir","bmi","waist","height","weight"]:
    df[c] = df.groupby("cycle")[c].transform(lambda s: s.fillna(s.median()))
for c in ["female","educ_hs","smoke_ever","smoke_now","race"]:
    df[c] = df.groupby("cycle")[c].transform(
        lambda s: s.fillna(s.mode().iloc[0] if not s.mode().empty else 0))


def nb(y, p, w, t):
    W = np.sum(w); flag = p >= t
    return (np.sum(w[flag & (y == 1)]) / W) - (np.sum(w[flag & (y == 0)]) / W) * t / (1 - t)


def fit(X, y, w, family):
    m = (LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=2000)
         if family == "lr" else
         HistGradientBoostingClassifier(max_iter=200, learning_rate=0.06,
                                        max_leaf_nodes=15, random_state=0))
    m.fit(X, y, sample_weight=w)
    return m


rows = []
for task in TASKS:
    base = df[df.cycle.isin(BASE) & df[task].notna()]
    prev_u = float(base[task].mean())
    if not (prev_u >= 0.01 and base[task].sum() >= 100):
        continue
    thr = max(round(prev_u, 2), 0.01)
    predsets = {"Q": Q} if task in Q_ONLY else {"Q": Q, "QE": QE}
    # B: hypertension is additionally followed into 2021-2023 here
    evals = list(EVAL)
    for ps_name, feats in predsets.items():
        Xb = base[feats].to_numpy(float); yb = base[task].to_numpy(float)
        wb = base.WTMEC2YR.to_numpy(float); wb = wb / wb.mean()
        for family in ["lr", "gbt"]:
            model = fit(Xb, yb, wb, family)
            for cyc in evals:
                ev = df[(df.cycle == cyc) & df[task].notna()]
                if len(ev) < 200:
                    continue
                ye = ev[task].to_numpy(float)
                we = ev.WTMEC2YR.to_numpy(float); we = we / we.mean()
                p0 = model.predict_proba(ev[feats].to_numpy(float))[:, 1]
                prev_d = float(np.sum(we * ye) / np.sum(we))
                rows.append(dict(
                    task=task, predset=ps_name, family=family, cycle=cyc,
                    midyear=MIDYEAR[cyc], prev_design=prev_d, thr=thr,
                    nb_fixed=nb(ye, p0, we, thr),
                    nb_half=nb(ye, p0, we, max(thr / 2, 0.005)),
                    nb_double=nb(ye, p0, we, min(thr * 2, 0.95)),
                    nb_tracking=nb(ye, p0, we, max(round(prev_d, 2), 0.01)),
                    thr_tracking=max(round(prev_d, 2), 0.01)))

m = pd.DataFrame(rows)
m.to_csv(f"{OUT}/threshold_mechanism.csv", index=False)
KEY = ["task", "predset", "family"]


def obsolete(d, col, k=2, ret=0.5, drop_2021_htn=True):
    """count combinations whose net benefit in `col` falls to <= ret x its initial value
    for k consecutive cycles, following the primary rule."""
    d = d.sort_values("midyear")
    if drop_2021_htn:
        d = d[~((d.task.isin(STOP_2017)) & (d.cycle == "2021-2023"))]
    n = tot = 0
    for key, s in d.groupby(KEY):
        tot += 1
        b = (s[col] <= ret * s[col].iloc[0]).to_numpy()
        if any(b[i] and b[i + 1] for i in range(len(b) - 1)):
            n += 1
    return n, tot


print("=" * 72)
print("A. PRESPECIFICATION v1 section 4 - threshold sensitivity (0.5x, 1x, 2x)")
print("=" * 72)
for col, lab in [("nb_half", "0.5x threshold"), ("nb_fixed", "1x (primary)"),
                 ("nb_double", "2x threshold")]:
    n, tot = obsolete(m, col)
    tasks = sorted({k[0] for k, s in m[~((m.task.isin(STOP_2017)) &
                                         (m.cycle == "2021-2023"))].sort_values("midyear").groupby(KEY)
                    for b in [(s[col] <= 0.5 * s[col].iloc[0]).to_numpy()]
                    if any(b[i] and b[i + 1] for i in range(len(b) - 1))})
    print(f"  {lab:16s} obsolete {n}/{tot} ({100*n/tot:.1f}%)   failing tasks: {tasks}")

print()
print("=" * 72)
print("B. PRESPECIFICATION v1 section 9 - hypertension with 2021-2023 included")
print("=" * 72)
h = m[m.task == "y_htn_measured"].sort_values("midyear")
n_stop, _ = obsolete(h, "nb_fixed", drop_2021_htn=True)
n_incl, _ = obsolete(h, "nb_fixed", drop_2021_htn=False)
print(f"  stopping at 2017-2018 (pre-specified): {n_stop}/4 combinations obsolete")
print(f"  including 2021-2023                  : {n_incl}/4 combinations obsolete")
print("  net benefit and design-based prevalence by cycle:")
for cyc, s in h.groupby("cycle", sort=False):
    print(f"    {cyc}  prev {s.prev_design.mean():.4f}  nb {s.nb_fixed.min():+.4f} to {s.nb_fixed.max():+.4f}")

print()
print("=" * 72)
print("C. Mechanism test - same frozen predictions, threshold tracks prevalence")
print("=" * 72)
for lab, col in [("threshold fixed at base period (as reported)", "nb_fixed"),
                 ("threshold tracks contemporaneous prevalence", "nb_tracking")]:
    n, tot = obsolete(m, col)
    print(f"  {lab:46s} obsolete {n}/{tot}")
h2 = h[h.cycle != "2021-2023"]
print("\n  measured hypertension, by cycle:")
print("    cycle      prev    thr_fixed  nb_fixed(min..max)   thr_track  nb_track(min..max)")
for cyc, s in h2.groupby("cycle", sort=False):
    print(f"    {cyc}  {s.prev_design.mean():.4f}   {s.thr.mean():.2f}      "
          f"{s.nb_fixed.min():+.4f}..{s.nb_fixed.max():+.4f}     {s.thr_tracking.mean():.2f}      "
          f"{s.nb_tracking.min():+.4f}..{s.nb_tracking.max():+.4f}")
r0 = h2[h2.cycle == "2003-2004"]
print(f"\n  fraction of initial net benefit retained at the 2009-2010 trough:")
tr = h2[h2.cycle == "2009-2010"]
for key in sorted(set(zip(h2.predset, h2.family))):
    a = r0[(r0.predset == key[0]) & (r0.family == key[1])]
    b = tr[(tr.predset == key[0]) & (tr.family == key[1])]
    print(f"    {key[0]:3s}/{key[1]:3s}  fixed {float(b.nb_fixed.iloc[0])/float(a.nb_fixed.iloc[0]):+.3f}   "
          f"tracking {float(b.nb_tracking.iloc[0])/float(a.nb_tracking.iloc[0]):+.3f}")
