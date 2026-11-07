# Statistical integrity report

Prepared 6 November 2026 for the manuscript *Calibration drift is not obsolescence: a
20-year decision-analytic study of 36 frozen clinical prediction models*.

This report states, in one place, every analytical decision that a reviewer could
reasonably contest, what was done about it, and where in the archive it can be checked.
It is written to be read alongside `PRESPECIFICATION_v1.md`,
`PRESPECIFICATION_v2_posthoc.md` and `RESULTS_MEMO_01.md`.

---

## 1. Pre-specification and its failure

| | |
|---|---|
| Plan written before any predictor–outcome association | Yes |
| SHA-256 of v1 plan | `5daf3197b26be2714c3f5f2b6fc240873792ce5f487360869f677e598af4aac8` |
| Work preceding the plan | file inventory; a recode-and-merge script printing row counts and missingness only |
| Falsification conditions pre-specified | Yes (§10 of v1) |
| Were they met | **Yes** |

The pre-specified acceptability boundary — net benefit at or below the better of treat-all
and treat-none — is near-vacuous when the decision threshold equals outcome prevalence,
because the net benefit of treat-all is then exactly zero by construction. The boundary
therefore collapsed to "worse than doing nothing", and almost nothing crossed it.

This was not concealed or reframed. The v1 result is reported in the main text
(`11_survival.py` reproduces it), a second specification was written and labelled post-hoc
throughout (`999a05d1a263f59d3150e40fa9cc3c49258b0f16e1d9cbeecdb39d522a9e1777`), and every
result in the manuscript carries a pre-specified or post-hoc label. A paper whose argument
is that acceptability boundaries are chosen carelessly is in no position to hide one.

**The primary acceptability boundary (50% retention of initial net benefit) is post-hoc.**
This is stated in the Abstract, the Methods and the Discussion.

---

## 2. The clustering problem — the single largest limitation

The four model–task combinations that lost clinical value are four variants of **one**
prediction task (measured hypertension: two predictor sets × two model families). At the
level of independent failure events the study observes **n = 1**, not n = 4.

Consequences, all of which are stated in the manuscript:

- The positive predictive values 0.17 (4/23) and 0.15 (3/20) are **descriptive**, not
  inferential, and **no interval is attached to either**. An earlier draft reported
  Clopper–Pearson intervals (0.05–0.39 and 0.03–0.38) and said they were wide because of
  the clustering. That was wrong in the paper's own favour: a Clopper–Pearson interval is a
  function of 4 and 23 alone and is blind to clustering, and clustering makes such an
  interval **too narrow** for the information present, not too wide. The intervals and the
  sentence have both been removed.
- No hypothesis test is performed on the signal–utility correspondence, and none should be.
- The task-level figures are also given (calibration-in-the-large fired on 8 of 10 tasks,
  slope on 9 of 10, against 1 task that lost value), so a reader can see the comparison at
  the level where the units are independent.
- The Abstract and the Conclusion both carry the caveat.

**This limits the paper to a demonstration that the discordance can be large, not an
estimate of how large it is in general.** The manuscript makes that claim and no stronger.

---

## 3. Estimator consistency: the decision threshold

Every evaluation quantity — AUROC, calibration-in-the-large, calibration slope, net
benefit, prevalence — is weighted by `WTMEC2YR`. Point estimates use the weights only;
`SDMVSTRA` and `SDMVPSU` enter the design-based confidence intervals computed by
`40_bootstrap_ci.py`, not the point estimates. The
decision threshold, however, is the **unweighted** base-period sample proportion rounded to
0.01, because that is the conventional default and is computable from base-period data
alone.

NHANES oversamples, so the two estimators differ. The threshold sits **above** the
design-based base-period prevalence for all ten tasks, by 0.007 to 0.070 (largest:
measured hypertension, 0.49 against 0.420, of which rounding explains only 0.0037).

The Methods state this explicitly and do **not** claim the threshold is the break-even
point for the target population.

**Sensitivity analysis** (`30_fit_weighted_threshold.py`, outputs in `out_wtd/`): the whole
follow-up repeated with design-based thresholds, nothing else changed.

| Quantity | Primary (unweighted threshold) | Sensitivity (weighted threshold) |
|---|---|---|
| Combinations triggered: CITL / slope / utility / AUROC | 23 / 20 / 4 / 0 | 23 / 20 / 4 / 0 |
| False alarms: CITL / slope | 19 / 17 | 19 / 17 |
| Positive predictive value: CITL / slope | 0.17 / 0.15 | 0.17 / 0.15 |
| Tasks losing clinical value | 1 (measured hypertension) | 1 (measured hypertension) |
| Apparent 20-year survival: utility / AUROC / CITL / slope | 89% / 100% / 36% / 44% | 89% / 100% / 36% / 44% |
| **Utility breach rate across evaluation points** | **6.6%** | **4.1%** |
| **Timing of the hypertension failure** | **6 years after freezing** | **8 years** |
| **Secondary boundary: triggers, true, PPV** | **3, 1, 0.33** | **2, 0, 0.00** |
| **Retention 0.25 / 0.75 cells** | **4/36 and 13/36** | **0/36 and 11/36** |

The central comparison — trigger counts, false alarms, positive predictive values, the
single failing task and the four 20-year survival figures — is identical. Four quantities
move, and the manuscript reports all four. **The signal–value discordance does not depend
on how the threshold was derived; the utility rule is itself threshold-sensitive at the
margins, which the retention-fraction analysis shows from a different direction.**

An earlier draft of this report and of the manuscript claimed that *every* reported
quantity was unchanged except the breach rate. That was wrong and was caught by an
independent audit; the four moving quantities above are the correction.

*Open decision for the corresponding author:* whether the design-based threshold should
become the primary analysis with the unweighted one as sensitivity, rather than the reverse.
The manuscript currently presents the unweighted threshold as primary, on the ground that it
is what a deploying team would actually compute. Either choice is defensible; the result is
the same.

---

## 4. The one genuine failure, and what was wrong with the first explanation

The measured-hypertension failure was initially attributed to the 2017 ACC/AHA
redefinition of hypertension. **That explanation was wrong and was removed.** Three facts
falsify it:

1. Net benefit bottoms at the 2009–2010 cycle; the guideline is from 2017.
2. The outcome was defined identically at ≥130/≥80 mmHg in every cycle, including before
   2017. No definition changed in this study.
3. Prevalence *fell* rather than rose.

The corrected account: the threshold was fixed at 0.49 from the base period (unweighted
0.4863 rounded up). Design-based prevalence then fell from 42.0% to 34.4% at 2009–2010 — an
18% relative decline — before recovering to 40.2% by 2017–2018. At t = 0.49 each false
positive costs 0.96 true positives, so a falling event rate removes net benefit quickly.
The model, its predictors and the outcome definition were all unchanged.

This also explains why full refitting did not rescue it, and why policy P4 produced one
*additional* failure (`y_mi`/QE/gbt): the loss is in the decision rule, not in the fitted
function.

This was originally an inference from the co-movement of prevalence and net benefit. It is
now a test. `31_threshold_mechanism.py` holds the frozen predictions exactly as they are and
lets only the decision threshold track contemporaneous weighted prevalence. Obsolescence
falls from 4 of 36 combinations to 1 of 36, and at the 2009–2010 trough the hypertension
models retain 72–82% of their initial net benefit instead of −17% to +11%. Nothing about the
model changed; only the threshold did.

Three references cited only in support of the discarded explanation (Whelton, Muntner,
Meyer-Eschenbach) were removed from the reference list.

---

## 5. Coding errors found and corrected before any result was reported

| # | Error | Effect | Fix |
|---|---|---|---|
| 1 | Intercept-only recalibration (P1) paired the intercept from a two-parameter fit with a slope fixed at 1.0 — not a valid intercept-only fit | P1 appeared to make models worse (14/36 obsolete) | Replaced with a weighted offset fit (`calib_intercept_only`). P1 now gives 4/36, identical to P0, P2 and P3. The erroneous figure was never reported. |
| 2 | The secondary utility boundary was coded as "net benefit falls by ≥1 net TP/1000 from its own initial value" instead of the specified "adds <1 net TP/1000 over the better default" | 36.9% instead of 5.0% | `20_tables.py` now implements the specification. |
| 3 | Figure 1 carried hardcoded sample sizes | Wrong n displayed | Sizes now computed from the data: 8,918 base / 48,034 evaluation / 56,952 total. |
| 4 | Figure 2 plotted relative net benefit on a shared scale, reaching 5× | Visually misleading | Switched to absolute net benefit with per-panel scales. |
| 5 | Figure 1 heading read "10 tasks × 2 predictor sets × 2 model families = 36" | Arithmetic is 40 | Corrected: the two obesity tasks take the questionnaire set only, since an outcome is never its own predictor. 8×4 + 2×2 = 36. |
| 6 | `SDMVSTRA` and `SDMVPSU` were read into local variables and never used; the Methods nevertheless described every quantity as "design-based", and the pre-specified bootstrap confidence intervals (§8) were never computed | A methods claim contradicted by the archived code | `40_bootstrap_ci.py` now performs the specified PSU-within-stratum bootstrap (500 replicates, seed 20261106) and asserts its point estimates reproduce the primary run to 1e-15. The Methods now say exactly what is weighted and what is design-based. |
| 7 | Two pre-specified sensitivity analyses (§4 thresholds at 0.5× and 2×; §9 hypertension including 2021–2023) were silently dropped | Undisclosed deviation from the plan | Both run and reported (`31_threshold_mechanism.py`): 0/36 and 11/36; all four variants obsolete either way. |
| 8 | Results said the threshold exceeded weighted prevalence "at the first evaluation cycle" for all ten tasks | True of the base period (10/10), false of the first evaluation cycle (6/10) | Corrected to "base period". |
| 9 | 212/320 printed as 66.2% | 66.25% rounds to 66.3% | Corrected. |
| 10 | The claim that no statistic computed from a model's own predictions can detect a threshold–prevalence gap | False: calibration-in-the-large is such a statistic, and it fired on these models | Replaced throughout with the sharper and true claim — that a band breach cannot distinguish a drifted population from a decayed model, and the two demand opposite responses. |
| 11 | `40_bootstrap_ci.py` derived its per-cell seed from `hash()` of a tuple of strings. Python salts `hash()` per process unless `PYTHONHASHSEED` is set, so the replicate draws were **not** reproducible across runs — contrary to what the pre-specification and the manuscript both claim | A reproducibility claim contradicted by the code | Seed is now a SHA-256 digest of the cell key. The whole bootstrap was re-run with the fixed seed. **Point estimates are unchanged** (49.4% and 45.3% breach rates), as the end-of-run assertion against `followup_metrics.csv` confirms. The interval-based percentages moved slightly and the manuscript now reports the new values: whole 95% interval outside the band at **33.4%** of points for slope (was 34.4%) and **25.9%** for calibration-in-the-large (was 25.6%); **32%** of slope breaches not separable from sampling error (was 30%); 43% for calibration-in-the-large, unchanged. |
| 12 | The exact 95% intervals attached to the two positive predictive values were justified by saying the four failures are variants of one task, so the intervals are wide | The justification ran backwards: dependence among the 36 units makes a binomial interval computed over them **too narrow**, not too wide | Both intervals removed. The manuscript states in the Results and again in the Limitations that the positive predictive values are descriptive and carry no interval. |
| 13 | `20_tables.py` bound `tot` to 320 evaluation points, then rebound it to 36 model–tasks inside a later loop, so `key_numbers.json` recorded `n_eval_points: 36` | A derived-data file disagreeing with the manuscript | The script is fixed, not just the JSON: a clean re-run now writes 320. Typographic normalisation of the tables was moved into the same script, and its output is byte-identical to the hand-corrected tables. |
| 14 | Figure 1 stated that net benefit, calibration-in-the-large, calibration slope and AUROC were "design-based (NHANES MEC weights, strata, PSUs)" | The same false claim as error 6, corrected in the Methods but left standing in the figure | Figure 1 regenerated from `21_figures.py`; it now reads "weighted by the NHANES MEC weight". Its alt text was corrected to match. Regenerating produced the other five figures byte-identical. |

Errors 1 and 2 are also recorded in the archive README. Every one of them was found by a
verification pass run against the code and the data, not by inspection of the prose.

---

## 6. Claims checked and withdrawn or softened

| Original claim | Status | Replacement |
|---|---|---|
| Hypertension failure caused by outcome-definition drift | **Withdrawn** — falsified (§4) | Fixed threshold, moving prevalence |
| "30% relative decline" in hypertension prevalence | **Corrected** — mixed an unweighted base with a design-based later cycle | about 18%, design-based throughout |
| Threshold is the break-even point | **Withdrawn** — it is not, for any task | Threshold sits above design-based prevalence by 0.007–0.070; sensitivity analysis provided |
| "Standardized net benefit" | **Corrected** — the formula is unstandardized | "net benefit" |
| Calibration slope "fired on 3 of 4 at the same time" | **Corrected** — misread a median | it fired on 3 of the 4, at 4, 6 and 14 years since the models were frozen |
| Poverty–income ratio "the only predictor missing in more than 5%" | **Corrected** — waist circumference is 5.71% | "the most incomplete predictor at 9.5%" |
| "Replication" of Williams et al. | **Softened** — that study is a simulation | "consistent with the direction of both" |
| Obsolescence has "no empirical literature" | **Scoped** | "the phase of that lifecycle with the thinnest empirical literature" |
| No framework monitors prevalence against the threshold | **Scoped** | "we are not aware of one that…" |
| "No study has examined" the signal–utility correspondence | **Scoped twice** | now "reported only incidentally, never across a portfolio over calendar time" — Drebin et al. 2025 reports it for one task family over 22 years |
| Obsolescence has "no empirical literature" | **Withdrawn** | "the least evidenced phase", with the two empirical studies that do exist cited |
| The apparent-lifespan-depends-on-the-rule idea presented as new | **Scoped** | the generic proposition is credited to Kim et al. 2025; only the magnitude for fitness for use is claimed as new |
| Calibration drifts while AUROC holds, presented as a finding | **Reframed as replication** | credited to Davis et al. 2017 and two 2026 replications; only the net-benefit half is claimed as new |
| The low rate of utility loss implicitly generalised | **Scoped** | stated in Limitations as a property of this portfolio, with Drebin et al. as the counterexample |

---

## 7. What the study cannot support, stated plainly

- **The mechanism claim is now tested, not inferred** (§4). What remains untested is whether
  other time-varying factors contributed; the test shows the threshold alone is sufficient to
  account for the loss, not that nothing else moved.
- **No causal claim about patients.** These models were never deployed. The study observes
  which signal a governance committee *could* have acted on, not what acting would have
  achieved. The Conclusion says this.
- **No generalization of the discordance magnitude**, for the clustering reason in §2.
- **No claim about acute, care-sensitive outcomes.** NHANES outcomes are slow-moving
  population phenomena. The Discussion states that an electronic health record cohort is
  the necessary replication, and predicts the discordance would be larger, not smaller —
  as a prediction, not a finding.
- **No claim that calibration monitoring is useless.** It detected the one genuine failure,
  and early. The claim is about its positive predictive value as a *retirement* trigger.
- **No inference from a single survey programme to all populations.** NHANES is United
  States, civilian, non-institutionalized, adults ≥20.

---

## 8. Reproducibility

Everything reported is reproducible from the archive at
`doi:10.5281/zenodo.23988406`, mirrored on GitHub. Scripts run in numeric order; the
numbering preserves the order in which they were written, so that everything below `10_`
is visibly structural and predates the pre-specification. `SHA256SUMS` covers every file;
`file_inventory.csv` carries the SHA-256 of each of the 77 NHANES `.xpt` source files, which
are not redistributed.

Raw NHANES data are public domain and downloadable from
https://wwwn.cdc.gov/nchs/nhanes/ .
