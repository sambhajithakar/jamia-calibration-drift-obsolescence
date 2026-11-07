"""
Manuscript tables. Recomputes every reported number from followup_metrics.csv
under the PRESPECIFICATION_v2_posthoc.md boundary, so the tables do not depend
on intermediate files.

Boundary (post-hoc, v2): a model is obsolete at the first evaluation cycle at
which NB_model <= RETENTION * NB_ref, sustained for k consecutive cycles, where
NB_ref is the model's net benefit at the first evaluation cycle (2003-2004).
"""
import numpy as np, pandas as pd

OUT = "/home/claude/jamia2026/03_analysis/out"
TAB = "/home/claude/jamia2026/04_manuscript/tables"
import os; os.makedirs(TAB, exist_ok=True)

BASE_END = 2001.5
RET, K = 0.50, 2
CITL_BAND  = (-0.2, 0.2)
SLOPE_BAND = (0.8, 1.25)
AUROC_DROP = 0.05

res = pd.read_csv(f"{OUT}/followup_metrics.csv")
KEY = ["task", "predset", "family"]

# ---- reference quantities from the first evaluation cycle, policy P0 ----
p0 = res[res.policy == "P0"].sort_values("midyear")
ref = p0.groupby(KEY).agg(nb_ref=("nb", "first"), auroc0=("auroc", "first")).reset_index()
res = res.merge(ref, on=KEY, how="left")

# ---- rule flags ----
res["fail_util"]  = res.nb <= RET * res.nb_ref
res["fail_abs"]   = res.nb <= res.nb_default + 0.001        # per PRESPECIFICATION_v2 secondary boundary
res["fail_citl"]  = ~res.citl.between(*CITL_BAND)
res["fail_slope"] = ~res.slope.between(*SLOPE_BAND)
res["fail_auroc"] = res.auroc <= (res.auroc0 - AUROC_DROP)

RULES = [("util", "Net benefit falls below 50% of its initial value"),
         ("abs",  "Net benefit adds <1 net true positive per 1000 over the default strategy"),
         ("citl", "Calibration-in-the-large outside [-0.20, 0.20]"),
         ("slope","Calibration slope outside [0.80, 1.25]"),
         ("auroc","AUROC falls >0.05 below its initial value")]


def first_sustained(g, col, k=K):
    g = g.sort_values("midyear")
    f, y = g[col].to_numpy(), g["midyear"].to_numpy()
    for i in range(len(f) - k + 1):
        if f[i:i + k].all():
            return y[i]
    return None


def event_table(df, col, k=K):
    """per model-task: event indicator and time since freeze"""
    rec = []
    for key, g in df.groupby(KEY):
        t = first_sustained(g, col, k)
        last = g.midyear.max()
        rec.append(dict(zip(KEY, key), event=int(t is not None),
                        dur=(t if t is not None else last) - BASE_END))
    return pd.DataFrame(rec)


# ======================= TABLE 1 — the portfolio =======================
elig = pd.read_csv(f"{OUT}/task_eligibility.csv")
ana = pd.read_csv(f"{OUT}/analytic.csv.gz", usecols=["cycle"])
n_mt = p0.groupby(KEY).ngroups

t1 = elig.rename(columns={"task": "Outcome", "base_n": "Base-period n",
                          "base_events": "Base-period events",
                          "base_prev": "Base-period prevalence",
                          "eligible": "Entered portfolio"})
# cycles followed, and model-task count, per task
foll = (p0.groupby("task")
          .agg(**{"Evaluation cycles": ("cycle", "nunique"),
                  "Last cycle": ("cycle", "max")}).reset_index())
nmt = p0.groupby("task").apply(lambda g: g.groupby(["predset", "family"]).ngroups,
                               include_groups=False).rename("Model-task combinations").reset_index()
t1 = (t1.merge(foll, left_on="Outcome", right_on="task", how="left")
        .merge(nmt, left_on="Outcome", right_on="task", how="left")
        .drop(columns=[c for c in ["task_x", "task_y", "task"] if c in t1.columns or True],
              errors="ignore"))
t1["Decision threshold"] = t1["Outcome"].map(
    p0.groupby("task").thr.first().to_dict())
LABEL_KEYS = list(t1["Outcome"])
LABEL = {"y_undiag_diabetes": "Undiagnosed diabetes", "y_obesity": "Obesity",
         "y_central_obesity": "Central obesity", "y_htn_measured": "Measured hypertension",
         "y_arthritis": "Arthritis", "y_chf": "Heart failure",
         "y_chd": "Coronary heart disease", "y_mi": "Myocardial infarction",
         "y_stroke": "Stroke", "y_cancer": "Cancer"}
t1["Outcome"] = t1["Outcome"].map(LABEL)
t1["Entered portfolio"] = t1["Entered portfolio"].map({True: "Yes", False: "No"})
# weighted base-period prevalence, shown alongside the unweighted one that sets the threshold
_an = pd.read_csv(f"{OUT}/analytic.csv.gz")
_base = _an[_an.cycle.isin(["1999-2000", "2001-2002"])]
_wprev = {}
for _tk in LABEL_KEYS:
    _b = _base[_base[_tk].notna()]
    _wprev[LABEL[_tk]] = float(np.average(_b[_tk].to_numpy(float),
                                          weights=_b.WTMEC2YR.to_numpy(float)))
t1["Base-period prevalence, weighted"] = (
    t1["Outcome"].map(_wprev) * 100).round(1).astype(str) + "%"
t1["Base-period prevalence"] = (t1["Base-period prevalence"] * 100).round(1).astype(str) + "%"
t1["Decision threshold"] = (t1["Decision threshold"] * 100).round(0).astype(int).astype(str) + "%"
t1 = t1[["Outcome", "Base-period n", "Base-period events", "Base-period prevalence",
         "Base-period prevalence, weighted",
         "Decision threshold", "Evaluation cycles", "Last cycle",
         "Model-task combinations", "Entered portfolio"]]
t1.to_csv(f"{TAB}/table1_portfolio.csv", index=False)

print("=" * 78)
print("TABLE 1 — prediction tasks and portfolio composition")
print("=" * 78)
print(t1.to_string(index=False))

# ============ TABLE 2 — how often each monitoring signal fires ============
rows = []
pp = res[res.policy == "P0"]
tot = len(pp)
for col, label in RULES:
    ev = event_table(pp, f"fail_{col}")
    rows.append(dict(Signal=label,
                     **{f"Evaluation points breached, n/{tot} (%)":
                            f"{int(pp[f'fail_{col}'].sum())} ({100 * pp[f'fail_{col}'].mean():.1f})",
                        f"Model-task combinations triggered, n/{len(ev)} (%)":
                            f"{int(ev.event.sum())} ({100 * ev.event.mean():.1f})"}))
t2 = pd.DataFrame(rows)
t2.to_csv(f"{TAB}/table2_breach_rates.csv", index=False)

print()
print("=" * 78)
print("TABLE 2 — frequency with which each monitoring signal fires (P0, frozen)")
print("=" * 78)
print(t2.to_string(index=False))

# ======== TABLE 3 — agreement with loss of clinical usefulness ========
ev_util = event_table(pp, "fail_util").set_index(KEY)
rows = []
for col, label in RULES:
    if col == "util":
        continue
    ev_r = event_table(pp, f"fail_{col}").set_index(KEY)
    j = ev_util.join(ev_r, lsuffix="_u", rsuffix="_r")
    lost  = j.event_u == 1          # genuinely lost half its value
    kept  = j.event_u == 0          # never did
    fired = j.event_r == 1
    lead = (j.loc[lost & fired, "dur_r"] - j.loc[lost & fired, "dur_u"])
    rows.append({
        "Signal": label,
        "Fires on models that lost half their value": f"{int((lost & fired).sum())}/{int(lost.sum())}",
        "False alarms (fires although value retained)": f"{int((kept & fired).sum())}/{int(kept.sum())}",
        "Median lead/lag vs utility loss, y":
            (("0.0" if abs(lead.median())<0.05 else f"{lead.median():+.1f}") if len(lead) else "\u2014"),
        "Positive predictive value for utility loss":
            (f"{(lost & fired).sum() / fired.sum():.2f}" if fired.sum() else "\u2014"),
    })
t3 = pd.DataFrame(rows)
t3.to_csv(f"{TAB}/table3_agreement.csv", index=False)

print()
print("=" * 78)
print("TABLE 3 — agreement between monitoring signals and loss of clinical usefulness")
print("=" * 78)
print(t3.to_string(index=False))

# ============== TABLE 4 — sensitivity and maintenance ==============
rows = []
for ret in [0.25, 0.50, 0.75]:
    tmp = pp.copy(); tmp["f"] = tmp.nb <= ret * tmp.nb_ref
    for k in [1, 2, 3]:
        ev = event_table(tmp, "f", k)
        rows.append(dict(Analysis=f"Retention {ret:.2f}, k={k}",
                         **{"Obsolete model-tasks": f"{int(ev.event.sum())}/{len(ev)}",
                            "%": round(100 * ev.event.mean(), 1)}))
for k in [1, 2, 3]:
    ev = event_table(pp, "fail_abs", k)
    rows.append(dict(Analysis=f"Secondary boundary (adds <1 net TP/1000), k={k}",
                     **{"Obsolete model-tasks": f"{int(ev.event.sum())}/{len(ev)}",
                        "%": round(100 * ev.event.mean(), 1)}))
for pol in ["P0", "P1", "P2", "P3", "P4"]:
    ev = event_table(res[res.policy == pol], "fail_util")
    rows.append(dict(Analysis=f"Maintenance policy {pol} (primary boundary)",
                     **{"Obsolete model-tasks": f"{int(ev.event.sum())}/{len(ev)}",
                        "%": round(100 * ev.event.mean(), 1)}))
# sensitivity: design-based decision thresholds (30_fit_weighted_threshold.py)
_wtd = "/home/claude/jamia2026/03_analysis/out_wtd/followup_metrics.csv"
if os.path.exists(_wtd):
    rw = pd.read_csv(_wtd)
    rw["fail_util"] = rw.nb <= 0.5 * rw.groupby(KEY).nb.transform("first")
    ev = event_table(rw[rw.policy == "P0"], "fail_util")
    rows.append(dict(Analysis="Weighted decision threshold (primary boundary)",
                     **{"Obsolete model-tasks": f"{int(ev.event.sum())}/{len(ev)}",
                        "%": round(100 * ev.event.mean(), 1)}))

# pre-specified threshold sensitivity (0.5x, 2x) and the mechanism test, from 31_
_mech = "/home/claude/jamia2026/03_analysis/out/threshold_mechanism.csv"
if os.path.exists(_mech):
    mm = pd.read_csv(_mech)
    mm = mm[~((mm.task == "y_htn_measured") & (mm.cycle == "2021-2023"))]
    for col, lab in [("nb_half", "Decision threshold halved (pre-specified)"),
                     ("nb_double", "Decision threshold doubled (pre-specified)"),
                     ("nb_tracking", "Threshold tracks contemporaneous prevalence")]:
        n = tot = 0
        for _k, s in mm.sort_values("midyear").groupby(KEY):
            tot += 1
            bb = (s[col] <= 0.5 * s[col].iloc[0]).to_numpy()
            if any(bb[i] and bb[i + 1] for i in range(len(bb) - 1)):
                n += 1
        rows.append(dict(Analysis=lab,
                         **{"Obsolete model-tasks": f"{n}/{tot}",
                            "%": round(100 * n / tot, 1)}))

t4 = pd.DataFrame(rows)
t4.to_csv(f"{TAB}/table4_sensitivity.csv", index=False)

print()
print("=" * 78)
print("TABLE 4 — sensitivity analyses and maintenance policies")
print("=" * 78)
print(t4.to_string(index=False))

# ---------- key numbers for the abstract, written to a file ----------
ev_u = event_table(pp, "fail_util")
key = {
    "n_model_tasks": int(n_mt),
    "n_eval_points": int(tot),
    "n_tasks": int(pp.task.nunique()),
    "pct_points_slope": round(100 * pp.fail_slope.mean(), 1),
    "pct_points_citl":  round(100 * pp.fail_citl.mean(), 1),
    "pct_points_util":  round(100 * pp.fail_util.mean(), 1),
    "pct_points_abs":   round(100 * pp.fail_abs.mean(), 1),
    "pct_points_auroc": round(100 * pp.fail_auroc.mean(), 1),
    "obsolete_util": int(ev_u.event.sum()), "n": len(ev_u),
    "tasks_lost_value": sorted(ev_u[ev_u.event == 1].task.unique().tolist()),
    "median_time_to_loss_y": float(ev_u[ev_u.event == 1].dur.median()) if ev_u.event.sum() else None,
}
pd.Series(key).to_json(f"{TAB}/key_numbers.json", indent=2)
print()
print("KEY NUMBERS")
for k, v in key.items():
    print(f"  {k:24s} {v}")
