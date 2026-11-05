"""
Useful life as a time-to-event outcome, determinants, maintenance-policy comparison,
and comparison of candidate retirement rules.

Executes PRESPECIFICATION_v1.md sections 5-7.
"""
import numpy as np, pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter

OUT = "/home/claude/jamia2026/03_analysis/out"
BASE_END = 2001.5          # last base-period midyear; clock starts here
K = 2                      # pre-specified consecutive cycles

res = pd.read_csv(f"{OUT}/followup_metrics.csv")

def first_sustained(g, flagcol, k=K):
    """first midyear at which flagcol is True for k consecutive evaluation cycles"""
    g = g.sort_values("midyear")
    f = g[flagcol].to_numpy()
    y = g["midyear"].to_numpy()
    for i in range(len(f) - k + 1):
        if f[i:i+k].all():
            return y[i]
    return None

# ---------- rule flags ----------
base_auc = (res[res.policy == "P0"].sort_values("midyear")
            .groupby(["task","predset","family"]).auroc.first().rename("auroc0"))
res = res.merge(base_auc, on=["task","predset","family"], how="left")
res["fail_nb"]    = res.nb <= res.nb_default
res["fail_auroc"] = res.auroc <= (res.auroc0 - 0.05)
res["fail_citl"]  = ~res.citl.between(-0.2, 0.2)
res["fail_slope"] = ~res.slope.between(0.8, 1.25)

# ---------- primary: useful life under P0 ----------
recs = []
for (task, ps, fam), g in res[res.policy == "P0"].groupby(["task","predset","family"]):
    last = g.midyear.max()
    t = first_sustained(g, "fail_nb")
    recs.append(dict(task=task, predset=ps, family=fam,
                     event=int(t is not None),
                     dur=(t if t is not None else last) - BASE_END,
                     base_prev=g.prev.iloc[0], auroc0=g.auroc0.iloc[0],
                     prev_change_decade=(g.sort_values("midyear").prev.iloc[-1]
                                         - g.sort_values("midyear").prev.iloc[0])
                                        / ((last - g.midyear.min()) / 10)))
S = pd.DataFrame(recs)
S.to_csv(f"{OUT}/useful_life_P0.csv", index=False)

km = KaplanMeierFitter().fit(S.dur, S.event)
med = km.median_survival_time_
print("=" * 64)
print("PRIMARY — useful life of frozen models (P0), portfolio of "
      f"{len(S)} model-task combinations")
print("=" * 64)
print(f"obsolete within the window : {S.event.sum()} / {len(S)}")
print(f"median useful life         : {med if np.isfinite(med) else 'not reached'} years")
try:
    ci = km.confidence_interval_survival_function_
    print(f"survival at 10 y           : {float(km.predict(10)):.2f}")
    print(f"survival at 20 y           : {float(km.predict(20)):.2f}")
except Exception:
    pass
print()
print("time-to-obsolescence by task (P0):")
print(S.groupby("task").agg(n=("dur","size"), events=("event","sum"),
                            median_dur=("dur","median")).to_string())

# ---------- determinants ----------
print()
print("=" * 64); print("DETERMINANTS (Cox)"); print("=" * 64)
X = S.copy()
X["predset_QE"] = (X.predset == "QE").astype(int)
X["family_gbt"] = (X.family == "gbt").astype(int)
cph_df = X[["dur","event","base_prev","auroc0","prev_change_decade",
            "predset_QE","family_gbt"]].dropna()
if cph_df.event.sum() >= 10:
    cph = CoxPHFitter(penalizer=0.1).fit(cph_df, "dur", "event")
    print(cph.summary[["coef","exp(coef)","exp(coef) lower 95%",
                       "exp(coef) upper 95%","p"]].round(3).to_string())
else:
    print(f"only {int(cph_df.event.sum())} events — Cox not estimated (pre-specified fallback)")

# ---------- maintenance policies ----------
print()
print("=" * 64); print("MAINTENANCE POLICIES"); print("=" * 64)
pol_rows = []
for pol in ["P0","P1","P2","P3","P4"]:
    rr = []
    for (task, ps, fam), g in res[res.policy == pol].groupby(["task","predset","family"]):
        last = g.midyear.max()
        t = first_sustained(g, "fail_nb")
        rr.append(dict(event=int(t is not None),
                       dur=(t if t is not None else last) - BASE_END))
    d = pd.DataFrame(rr)
    k2 = KaplanMeierFitter().fit(d.dur, d.event)
    m = k2.median_survival_time_
    pol_rows.append(dict(policy=pol, n=len(d), obsolete=int(d.event.sum()),
                         pct_obsolete=round(100*d.event.mean(),1),
                         median_life=(m if np.isfinite(m) else np.nan)))
P = pd.DataFrame(pol_rows); P.to_csv(f"{OUT}/policy_comparison.csv", index=False)
print(P.to_string(index=False))

# unrescuable fraction: fails under every policy
f_by = (res.groupby(["task","predset","family","policy"])
        .apply(lambda g: first_sustained(g, "fail_nb") is not None, include_groups=False)
        .unstack("policy"))
unres = f_by.all(axis=1)
print()
print(f"unrescuable (obsolete under ALL five policies): {int(unres.sum())} / {len(unres)}"
      f"  ({100*unres.mean():.1f}%)")
f_by.assign(unrescuable=unres).to_csv(f"{OUT}/unrescuable.csv")

# ---------- retirement rule comparison ----------
print()
print("=" * 64); print("RETIREMENT RULES — when would each fire?"); print("=" * 64)
cmp_rows = []
for (task, ps, fam), g in res[res.policy == "P0"].groupby(["task","predset","family"]):
    r = dict(task=task, predset=ps, family=fam)
    for rule, col in [("net_benefit","fail_nb"), ("auroc_drop","fail_auroc"),
                      ("citl","fail_citl"), ("slope","fail_slope")]:
        t = first_sustained(g, col)
        r[rule] = (t - BASE_END) if t is not None else np.nan
    cmp_rows.append(r)
C = pd.DataFrame(cmp_rows); C.to_csv(f"{OUT}/rule_comparison.csv", index=False)
print(C[["net_benefit","auroc_drop","citl","slope"]]
      .agg(["count","median"]).round(2).to_string())
print()
both = C.dropna(subset=["net_benefit"])
for rule in ["auroc_drop","citl","slope"]:
    sub = both.dropna(subset=[rule])
    if len(sub):
        lead = (sub[rule] - sub["net_benefit"])
        print(f"{rule:11s}: fires on {len(sub):2d}/{len(both)} of the models that lost "
              f"net benefit; median lag vs net-benefit rule = {lead.median():+.1f} years")
    else:
        print(f"{rule:11s}: never fires on any model that lost net benefit")
