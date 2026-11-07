"""
Manuscript figures (6, the JAMIA maximum). All numbers recomputed from
followup_metrics.csv under the v2 post-hoc boundary.
"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyArrowPatch

OUT = "/home/claude/jamia2026/03_analysis/out"
FIG = "/home/claude/jamia2026/04_manuscript/figures"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 8.5, "axes.titleweight": "bold",
    "axes.labelsize": 8, "legend.fontsize": 7.2,
    "xtick.labelsize": 7.2, "ytick.labelsize": 7.2,
    "figure.dpi": 400, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
})
# Okabe-Ito, colour-vision-deficiency safe
C_UTIL, C_CITL, C_SLOPE, C_AUROC = "#0072B2", "#D55E00", "#CC79A7", "#009E73"
C_GREY, C_LIGHT = "#4D4D4D", "#BFBFBF"

BASE_END, RET, K = 2001.5, 0.50, 2
KEY = ["task", "predset", "family"]
LABEL = {"y_undiag_diabetes": "Undiagnosed diabetes", "y_obesity": "Obesity",
         "y_central_obesity": "Central obesity", "y_htn_measured": "Measured hypertension",
         "y_arthritis": "Arthritis", "y_chf": "Heart failure",
         "y_chd": "Coronary heart disease", "y_mi": "Myocardial infarction",
         "y_stroke": "Stroke", "y_cancer": "Cancer"}

res = pd.read_csv(f"{OUT}/followup_metrics.csv")
p0g = res[res.policy == "P0"].sort_values("midyear")
ref = p0g.groupby(KEY).agg(nb_ref=("nb", "first"), auroc0=("auroc", "first")).reset_index()
res = res.merge(ref, on=KEY, how="left")
res["fail_util"]  = res.nb <= RET * res.nb_ref
res["fail_citl"]  = ~res.citl.between(-0.2, 0.2)
res["fail_slope"] = ~res.slope.between(0.8, 1.25)
res["fail_auroc"] = res.auroc <= (res.auroc0 - 0.05)
res["nb_rel"]     = res.nb / res.nb_ref
pp = res[res.policy == "P0"].copy()


def first_sustained(g, col, k=K):
    g = g.sort_values("midyear")
    f, y = g[col].to_numpy(), g["midyear"].to_numpy()
    for i in range(len(f) - k + 1):
        if f[i:i + k].all():
            return y[i]
    return None


# ============================ FIGURE 1 — design ============================
fig, ax = plt.subplots(figsize=(7.1, 2.5))
ax.set_xlim(1998.0, 2024.2); ax.set_ylim(0, 1); ax.axis("off")
cycles = [1999.5, 2001.5, 2003.5, 2005.5, 2007.5, 2009.5, 2011.5, 2013.5,
          2015.5, 2017.5, 2022.0]
names = ["1999-\n2000", "2001-\n2002", "2003-\n2004", "2005-\n2006", "2007-\n2008",
         "2009-\n2010", "2011-\n2012", "2013-\n2014", "2015-\n2016", "2017-\n2018",
         "2021-\n2023"]
ana_c = pd.read_csv(f"{OUT}/analytic.csv.gz", usecols=["cycle"]).cycle
N_BASE = int(ana_c.isin(["1999-2000", "2001-2002"]).sum())
N_EVAL = int(len(ana_c) - N_BASE)
ax.add_patch(Rectangle((1998.25, 0.48), 4.15, 0.26, fc="#D9E8F5", ec=C_UTIL, lw=0.9))
ax.text(2000.33, 0.61, f"Training\n(base period)\nn = {N_BASE:,} adults",
        ha="center", va="center", fontsize=5.6, color=C_UTIL, fontweight="bold",
        linespacing=1.3)
ax.add_patch(Rectangle((2002.6, 0.48), 21.0, 0.26, fc="#F5F5F5", ec=C_GREY, lw=0.9))
ax.text(2013.0, 0.61, "Evaluation period \u2014 models frozen, followed forward 20.5 years\n"
        f"9 cycles, n = {N_EVAL:,} adults", ha="center", va="center",
        fontsize=7.2, color=C_GREY, linespacing=1.45)
for x, nm in zip(cycles, names):
    inbase = x <= 2001.5
    ax.plot([x, x], [0.40, 0.48], color=C_UTIL if inbase else C_GREY, lw=0.9)
    ax.plot(x, 0.38, "o", ms=3.4, color=C_UTIL if inbase else C_GREY)
    ax.text(x, 0.30, nm, ha="center", va="top", fontsize=5.9,
            color=C_UTIL if inbase else C_GREY)
ax.annotate("", xy=(2019.6, 0.38), xytext=(2020.6, 0.38),
            arrowprops=dict(arrowstyle="-", color=C_LIGHT, lw=0.8, ls=(0, (1, 1))))
ax.text(2020.1, 0.20, "2019-2020\nnot released\n(COVID-19)", ha="center", va="top",
        fontsize=5.5, color=C_LIGHT, style="italic")
ax.plot([1998.6, 2023.6], [0.38, 0.38], color=C_LIGHT, lw=0.6, zorder=0)
ax.text(1998.6, 0.90, "10 prediction tasks  ×  2 model families  ×  up to 2 predictor sets  =  "
        "36 model–task combinations", fontsize=7.4, fontweight="bold", color="black")
ax.text(1998.6, 0.845, "(the two obesity tasks take the questionnaire set only, since an outcome "
        "is never its own predictor)", fontsize=6.2, style="italic", color=C_GREY)
ax.text(1998.6, 0.80, "At every evaluation cycle: net benefit, calibration-in-the-large, "
        "calibration slope, AUROC — weighted by the NHANES MEC weight",
        fontsize=6.8, color=C_GREY)
fig.savefig(f"{FIG}/figure1_design.png"); plt.close(fig)

# ============== FIGURE 2 — net benefit trajectories, per task ==============
SHORT = {"y_undiag_diabetes": "Undiagnosed\ndiabetes", "y_obesity": "Obesity",
         "y_central_obesity": "Central obesity", "y_htn_measured": "Hypertension\n(measured)",
         "y_arthritis": "Arthritis", "y_chf": "Heart failure",
         "y_chd": "Coronary heart\ndisease", "y_mi": "Myocardial\ninfarction",
         "y_stroke": "Stroke", "y_cancer": "Cancer"}
tasks = [t for t in LABEL if t in pp.task.unique()]
fig, axes = plt.subplots(2, 5, figsize=(7.1, 3.7), sharex=True)
for ax, task in zip(axes.ravel(), tasks):
    g = pp[pp.task == task]
    for _, gg in g.groupby(["predset", "family"]):
        gg = gg.sort_values("midyear")
        ax.plot(gg.midyear, gg.nb, color=C_UTIL, lw=0.9, alpha=0.8)
    ref_lo = (g.groupby(["predset", "family"]).nb_ref.first() * RET).min()
    ref_hi = (g.groupby(["predset", "family"]).nb_ref.first() * RET).max()
    ax.axhspan(ref_lo, ref_hi, color=C_UTIL, alpha=0.18, zorder=1)
    ax.axhline(ref_lo, color=C_UTIL, lw=0.7, ls=(0, (3, 2)), zorder=2)
    lo = min(0, g.nb.min()); hi = max(g.nb.max(), ref_hi)
    pad = 0.12 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    ax.axhline(0, color=C_LIGHT, lw=0.6, zorder=0)
    ax.set_title(SHORT[task], fontsize=6.6, linespacing=1.15)
    ax.set_xlim(2002.5, 2023); ax.set_xticks([2005, 2012, 2020])
    ax.set_xticklabels(["'05", "'12", "'20"])
    ax.tick_params(axis="y", labelsize=6.0)
for ax in axes[:, 0]:
    ax.set_ylabel("Net benefit", fontsize=7.2)
for ax in axes[1, :]:
    ax.set_xlabel("Calendar year", fontsize=7.2)
fig.legend(handles=[Line2D([], [], color=C_UTIL, lw=1.1,
                           label="one frozen model (4 per task, 2 for obesity measures)"),
                    Line2D([], [], color=C_UTIL, lw=0.8, ls=(0, (3, 2)),
                           label="50% of the model's initial (2003-2004) net benefit")],
           loc="lower center", bbox_to_anchor=(0.5, -0.10), ncol=2, frameon=False,
           fontsize=6.6)
fig.subplots_adjust(hspace=0.45, wspace=0.42)
fig.savefig(f"{FIG}/figure2_net_benefit.png"); plt.close(fig)

# ============== FIGURE 3 — calibration drifts while utility holds ==============
fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.3))
for ax, (col, lab, band, colr) in zip(
        axes, [("slope", "Calibration slope", (0.8, 1.25), C_SLOPE),
               ("citl", "Calibration-in-the-large", (-0.2, 0.2), C_CITL),
               ("nb_rel", "Net benefit (relative)", (RET, 10), C_UTIL)]):
    for _, gg in pp.groupby(KEY):
        gg = gg.sort_values("midyear")
        ax.plot(gg.midyear, gg[col], color=colr, lw=0.6, alpha=0.4)
    med = pp.groupby("midyear")[col].median()
    ax.plot(med.index, med.values, color="black", lw=1.6, label="portfolio median")
    ax.axhspan(band[0], band[1], color=colr, alpha=0.11, zorder=0,
               label="acceptable range")
    ax.set_title(lab); ax.set_xlabel("Calendar year")
    ax.set_xlim(2002.5, 2023)
axes[0].set_ylim(0, 2.2); axes[1].set_ylim(-1.3, 1.3); axes[2].set_ylim(-0.1, 1.8)
axes[0].set_ylabel("Value at evaluation cycle")
axes[2].legend(loc="lower left", frameon=False, fontsize=6.3)
for ax, pct in zip(axes, [49.4, 45.3, 6.6]):
    ax.text(0.97, 0.96, f"{pct}% of\nevaluation\npoints outside", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.3, color=C_GREY,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))
fig.savefig(f"{FIG}/figure3_calibration_vs_utility.png"); plt.close(fig)

# ================== FIGURE 4 — breach rates and agreement ==================
fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.9),
                         gridspec_kw=dict(width_ratios=[1, 1.3]))
ax = axes[0]
labs = ["Calibration\nslope", "Calibration-\nin-the-\nlarge", "Net\nbenefit\n(50% loss)", "AUROC\n(>0.05\ndrop)"]
vals = [pp.fail_slope.mean() * 100, pp.fail_citl.mean() * 100,
        pp.fail_util.mean() * 100, pp.fail_auroc.mean() * 100]
cols = [C_SLOPE, C_CITL, C_UTIL, C_AUROC]
b = ax.bar(labs, vals, color=cols, width=0.58)
ax.tick_params(axis="x", labelsize=5.6)
for lbl in ax.get_xticklabels():
    lbl.set_linespacing(1.1)
for r, v in zip(b, vals):
    ax.text(r.get_x() + r.get_width() / 2, v + 1.4, f"{v:.1f}%", ha="center", fontsize=7.4,
            fontweight="bold")
ax.set_ylabel("% of 320 evaluation points breached")
ax.set_ylim(0, 60); ax.set_title("How often each signal fires")

ax = axes[1]
ev = {}
for col in ["util", "citl", "slope", "auroc"]:
    rec = []
    for key, g in pp.groupby(KEY):
        rec.append(dict(zip(KEY, key), ev=int(first_sustained(g, f"fail_{col}") is not None)))
    ev[col] = pd.DataFrame(rec).set_index(KEY).ev
E = pd.DataFrame(ev)
lost, kept = E.util == 1, E.util == 0
rows = [("Calibration slope", C_SLOPE, "slope"), ("Calibration-in-the-large", C_CITL, "citl"),
        ("AUROC (>0.05 drop)", C_AUROC, "auroc")]
y = np.arange(len(rows))
tp = [int((lost & (E[c] == 1)).sum()) for _, _, c in rows]
fp = [int((kept & (E[c] == 1)).sum()) for _, _, c in rows]
fn = [int((lost & (E[c] == 0)).sum()) for _, _, c in rows]
ax.barh(y, tp, color="#1F6F43", height=0.5, label=f"correctly flags a model that lost value (n={int(lost.sum())})")
ax.barh(y, fp, left=tp, color="#E8A33D", height=0.5, label=f"false alarm — value retained (n={int(kept.sum())})")
for i, (t, f, n) in enumerate(zip(tp, fp, fn)):
    ax.text(t + f + 0.4, i, f"{t} true, {f} false   (missed {n})", va="center", fontsize=6.9)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows])
ax.set_xlabel("Model–task combinations triggered (of 36)")
ax.set_xlim(0, 41); ax.invert_yaxis()
ax.set_title("What the signal means when it fires")
ax.legend(loc="upper center", bbox_to_anchor=(0.42, -0.26), ncol=1,
          frameon=False, fontsize=6.3)
fig.subplots_adjust(wspace=0.62)
fig.savefig(f"{FIG}/figure4_breach_and_agreement.png"); plt.close(fig)

# ==================== FIGURE 5 — when each rule fires ====================
recs = []
for key, g in pp.groupby(KEY):
    r = dict(zip(KEY, key))
    for col in ["util", "citl", "slope", "auroc"]:
        t = first_sustained(g, f"fail_{col}")
        r[col] = (t - BASE_END) if t is not None else np.nan
    r["last"] = g.midyear.max() - BASE_END
    recs.append(r)
T = pd.DataFrame(recs)
T["lab"] = (T.task.map(LABEL) + "  " + T.predset + "/" + T.family.str.upper())
T = T.sort_values(["util", "citl", "lab"], na_position="last").reset_index(drop=True)
fig, ax = plt.subplots(figsize=(7.1, 5.2))
for i, r in T.iterrows():
    ax.plot([0, r["last"]], [i, i], color="#E8E8E8", lw=3.0, solid_capstyle="butt", zorder=0)
    for col, colr, mk, ms in [("citl", C_CITL, "o", 4.4), ("slope", C_SLOPE, "s", 3.9),
                              ("auroc", C_AUROC, "^", 4.4), ("util", C_UTIL, "D", 4.6)]:
        if not np.isnan(r[col]):
            ax.plot(r[col], i, mk, color=colr, ms=ms, mew=0, zorder=3)
ax.set_yticks(range(len(T))); ax.set_yticklabels(T.lab, fontsize=5.8)
ax.invert_yaxis()
ax.set_xlabel("Years since the model was frozen")
ax.set_xlim(-0.6, 21.5); ax.set_xticks([0, 4, 8, 12, 16, 20])
ax.set_title("First sustained trigger of each monitoring rule, by model–task combination")
ax.legend(handles=[
    Line2D([], [], ls="", marker="o", color=C_CITL, ms=4.4, label="calibration-in-the-large"),
    Line2D([], [], ls="", marker="s", color=C_SLOPE, ms=3.9, label="calibration slope"),
    Line2D([], [], ls="", marker="D", color=C_UTIL, ms=4.6, label="net benefit lost half its value"),
    Line2D([], [], ls="", marker="^", color=C_AUROC, ms=4.4, label="AUROC drop >0.05 (never fires)"),
    Line2D([], [], color="#E8E8E8", lw=3, label="follow-up observed"),
], loc="lower right", frameon=False, fontsize=6.6, ncol=1)
ax.grid(axis="x", color="#F0F0F0", lw=0.6); ax.set_axisbelow(True)
fig.savefig(f"{FIG}/figure5_trigger_timing.png"); plt.close(fig)

# ============== FIGURE 6 — useful life under each definition ==============
from lifelines import KaplanMeierFitter
fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.6), sharey=True)
ax = axes[0]
for col, colr, lab, ls in [("util", C_UTIL, "Net benefit lost half its value", "-"),
                           ("citl", C_CITL, "Calibration-in-the-large out of range", "-"),
                           ("slope", C_SLOPE, "Calibration slope out of range", "-"),
                           ("auroc", C_AUROC, "AUROC drop >0.05", (0, (4, 2)))]:
    d = T.copy()
    d["event"] = d[col].notna().astype(int)
    d["dur"] = d[col].fillna(d["last"])
    km = KaplanMeierFitter().fit(d.dur, d.event)
    km.plot_survival_function(ax=ax, color=colr, lw=1.5, ls=ls, ci_show=False, label=lab)
ax.set_xlabel("Years since freezing"); ax.set_ylabel("Proportion still in service")
ax.set_ylim(0, 1.02); ax.set_xlim(0, 21)
ax.set_title("Retirement rule determines apparent lifespan")
ax.legend(frameon=False, fontsize=6.3, loc="lower left")

ax = axes[1]
for ret, colr, lab in [(0.25, "#9ECAE1", "25% retention"), (0.50, C_UTIL, "50% retention (primary)"),
                       (0.75, "#08306B", "75% retention")]:
    tmp = pp.copy(); tmp["f"] = tmp.nb <= ret * tmp.nb_ref
    rec = []
    for key, g in tmp.groupby(KEY):
        t = first_sustained(g, "f")
        rec.append(dict(event=int(t is not None),
                        dur=(t if t is not None else g.midyear.max()) - BASE_END))
    d = pd.DataFrame(rec)
    km = KaplanMeierFitter().fit(d.dur, d.event)
    km.plot_survival_function(ax=ax, color=colr, lw=1.5, ci_show=False, label=lab)
ax.set_xlabel("Years since freezing"); ax.set_ylim(0, 1.02); ax.set_xlim(0, 21)
ax.set_title("Sensitivity to the utility-loss boundary")
ax.legend(frameon=False, fontsize=6.3, loc="lower left")
fig.savefig(f"{FIG}/figure6_useful_life.png"); plt.close(fig)

T.to_csv(f"{FIG}/../tables/trigger_timing_data.csv", index=False)
print("figures written:")
for f in sorted(os.listdir(FIG)):
    print(f"  {f}  {os.path.getsize(os.path.join(FIG,f))/1024:.0f} KB")
print()
print("base-period n:", "see table 1")
print(f"evaluation-period rows: {int(pp.n.sum() / 1)}  (sum of per-cycle n across tasks)")
