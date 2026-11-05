"""
v2 (POST-HOC) boundary: a model is obsolete when it has lost half the clinical
value it delivered at entry into service. Re-derives event times from the metrics
already computed under v1; no model is refitted.
"""
import numpy as np, pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter

OUT = "/home/claude/jamia2026/03_analysis/out"
BASE_END, K, RETAIN = 2001.5, 2, 0.50
res = pd.read_csv(f"{OUT}/followup_metrics.csv")

ref = (res[res.policy == "P0"].sort_values("midyear")
       .groupby(["task","predset","family"]).agg(nb_ref=("nb","first"),
                                                 auroc0=("auroc","first")).reset_index())
res = res.merge(ref, on=["task","predset","family"], how="left")
res["fail_v2"]    = res.nb <= RETAIN * res.nb_ref
res["fail_abs"]   = res.nb <= res.nb_default + 0.001
res["fail_citl"]  = ~res.citl.between(-0.2, 0.2)
res["fail_slope"] = ~res.slope.between(0.8, 1.25)
res["fail_auroc"] = res.auroc <= (res.auroc0 - 0.05)

def first_sustained(g, col, k=K):
    g = g.sort_values("midyear"); f = g[col].to_numpy(); y = g.midyear.to_numpy()
    for i in range(len(f)-k+1):
        if f[i:i+k].all(): return y[i]
    return None

def survtab(pol, col, k=K):
    rr = []
    for key, g in res[res.policy == pol].groupby(["task","predset","family"]):
        t = first_sustained(g, col, k); last = g.midyear.max()
        rr.append(dict(task=key[0], predset=key[1], family=key[2],
                       event=int(t is not None),
                       dur=(t if t is not None else last) - BASE_END,
                       nb_ref=g.nb_ref.iloc[0], auroc0=g.auroc0.iloc[0],
                       base_prev=g.prev.iloc[0]))
    return pd.DataFrame(rr)

S = survtab("P0", "fail_v2"); S.to_csv(f"{OUT}/useful_life_v2_P0.csv", index=False)
km = KaplanMeierFitter().fit(S.dur, S.event); med = km.median_survival_time_
print("="*66); print("v2 (POST-HOC) — useful life, frozen models, lose-half-the-value boundary")
print("="*66)
print(f"obsolete within window : {S.event.sum()} / {len(S)}  ({100*S.event.mean():.0f}%)")
print(f"median useful life     : {med if np.isfinite(med) else 'not reached'} years")
for h in (6, 10, 14, 20):
    print(f"  still fit at {h:2d} y     : {float(km.predict(h)):.2f}")
print()
print(S.groupby("task").agg(n=("dur","size"), events=("event","sum"),
                            median=("dur","median")).to_string())

print()
print("="*66); print("SENSITIVITY"); print("="*66)
for rt in (0.25, 0.50, 0.75):
    res["f"] = res.nb <= rt*res.nb_ref
    for kk in (1,2,3):
        d = survtab("P0","f",kk); m = KaplanMeierFitter().fit(d.dur,d.event).median_survival_time_
        print(f"  retain {rt:.2f}, k={kk}: obsolete {int(d.event.sum()):2d}/{len(d)}, "
              f"median {m if np.isfinite(m) else 'NR'}")

print()
print("="*66); print("MAINTENANCE POLICIES (v2 boundary)"); print("="*66)
rows=[]
for pol in ["P0","P1","P2","P3","P4"]:
    d = survtab(pol,"fail_v2"); m = KaplanMeierFitter().fit(d.dur,d.event).median_survival_time_
    rows.append(dict(policy=pol, obsolete=int(d.event.sum()), n=len(d),
                     pct=round(100*d.event.mean(),1),
                     median_life=(m if np.isfinite(m) else np.nan)))
P=pd.DataFrame(rows); P.to_csv(f"{OUT}/policy_comparison_v2.csv", index=False)
print(P.to_string(index=False))
fb = (res.groupby(["task","predset","family","policy"])
        .apply(lambda g: first_sustained(g,"fail_v2") is not None, include_groups=False)
        .unstack("policy"))
print(f"\nunrescuable under all five policies: {int(fb.all(axis=1).sum())}/{len(fb)}"
      f"  ({100*fb.all(axis=1).mean():.1f}%)")

print()
print("="*66); print("THE CENTRAL COMPARISON — breach rate per evaluation point (P0)")
print("="*66)
p0 = res[res.policy=="P0"]
for lab,c in [("calibration-in-the-large","fail_citl"),("calibration slope","fail_slope"),
              ("AUROC drop > 0.05","fail_auroc"),
              ("lost half its value (v2)","fail_v2"),
              ("adds <1 per 1000 (abs)","fail_abs")]:
    print(f"  {lab:28s}: {100*p0[c].mean():5.1f}%")

print()
print("rule firing times vs the v2 utility rule (model-task level, P0):")
C=[]
for key,g in p0.groupby(["task","predset","family"]):
    r=dict(task=key[0])
    for lab,c in [("v2","fail_v2"),("citl","fail_citl"),("slope","fail_slope"),("auroc","fail_auroc")]:
        t=first_sustained(g,c); r[lab]=(t-BASE_END) if t is not None else np.nan
    C.append(r)
C=pd.DataFrame(C); C.to_csv(f"{OUT}/rule_comparison_v2.csv", index=False)
print(C[["v2","citl","slope","auroc"]].agg(["count","median"]).round(1).to_string())
fired_v2 = C.v2.notna()
print(f"\nmodels where the utility rule ever fires : {fired_v2.sum()}/{len(C)}")
for c in ["citl","slope"]:
    早 = C[~fired_v2][c].notna().sum()
    print(f"  {c:5s} fires on {早}/{(~fired_v2).sum()} models that NEVER lost half their value"
          f"  -> false alarms")
