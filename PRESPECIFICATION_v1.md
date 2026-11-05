# Pre-specification v1
**Study:** Time to obsolescence of frozen clinical prediction models across a portfolio of tasks, NHANES 1999–2023
**Written:** 6 November 2026, before any predictor–outcome association was computed.

## Blindness statement

Everything done before this document was structural: a file and variable inventory (`01_harmonisation_audit.py`) and a recode-and-merge build (`02_build_analytic.py`). The build printed row counts, outcome non-missing counts and predictor missingness. **No model has been fitted, no outcome prevalence has been examined, and no predictor–outcome association has been computed.** The SHA-256 of this file is recorded in `PRESPECIFICATION_v1.sha256` and all analysis scripts are numbered from `10_` onward so that the ordering is visible in the repository.

## 1. Data

Eleven NHANES cycles: 1999-2000, 2001-2002, 2003-2004, 2005-2006, 2007-2008, 2009-2010, 2011-2012, 2013-2014, 2015-2016, 2017-2018, 2021-2023. The 2019-2020 cycle does not exist as a standalone public release (COVID-19) and is treated as a gap in the calendar axis, not closed.

Analytic sample: adults aged ≥ 20, examined (`RIDSTATR` = 2) with a positive MEC weight, not pregnant. n = 56,952.

## 2. Portfolio of prediction tasks

Candidate outcomes, all defined identically in every cycle:

| Outcome | Definition |
|---|---|
| `y_undiag_diabetes` | HbA1c ≥ 6.5% among adults answering "no" to `DIQ010` |
| `y_obesity` | BMI ≥ 30 |
| `y_central_obesity` | waist > 88 cm (women) or > 102 cm (men) |
| `y_htn_measured` | mean SBP ≥ 130 or mean DBP ≥ 80 |
| `y_arthritis`, `y_chf`, `y_chd`, `y_mi`, `y_stroke`, `y_cancer` | self-reported physician diagnosis (`MCQ160A/B/C/E/F`, `MCQ220`) |

**Inclusion rule, fixed now:** a task enters the portfolio if, in the base period, it has prevalence ≥ 1% **and** ≥ 100 events. Tasks failing this are reported in a table and excluded, because net benefit at a low threshold is unstable with few events.

**Predictor sets**
- **Q (questionnaire only):** age, sex, race/ethnicity, education ≥ high school, poverty–income ratio, ever smoked, currently smokes.
- **Q+E (adds examination):** + BMI, waist circumference, height, weight.

An outcome is never used as its own predictor; for `y_obesity` and `y_central_obesity` the Q+E set omits BMI, waist, height and weight, so those two tasks appear with the Q set only.

**Model families:** L2-penalised logistic regression, and gradient-boosted trees. Hyperparameters tuned once, in the base period only, by 5-fold cross-validation optimising log loss (not AUROC), and then frozen.

The portfolio is tasks × eligible predictor sets × model families. Target ≥ 30 model–task combinations. If fewer than 30 are eligible, the Cox model in §7 is dropped and useful life is reported descriptively; this fallback is fixed now, not after seeing the count.

## 3. Base period and freezing

Models are trained on the first two cycles (1999-2000 and 2001-2002) and then **frozen**. Every later cycle is an evaluation period. The base period is never re-entered.

## 4. Decision threshold

For each task the primary decision threshold is the **base-period prevalence of that outcome, rounded to the nearest 0.01**. This is the threshold at which treat-all and treat-none break even, it is the conventional default when no external utility is available, and it uses only base-period data, so it is blind to every evaluation cycle.

Sensitivity analyses at 0.5× and 2× the primary threshold.

## 5. Failure boundary — the definition of obsolescence

At each evaluation cycle compute the model's net benefit at the decision threshold *t*:

    NB_model = TP/N − (FP/N) × t/(1−t)

and the net benefit of the better default strategy, max(NB_treat-all, 0), where NB_treat-none = 0 by construction.

**A model is obsolete at the first evaluation cycle at which NB_model ≤ max(NB_treat-all, 0), sustained for k = 2 consecutive evaluation cycles.** The event time is the first of those cycles. Sensitivity analyses at k = 1 and k = 3.

Net benefit is used rather than AUROC or calibration because fitness for purpose is a decision-analytic property. A model can be well calibrated and clinically worthless, and whether that happens is the question the study exists to answer.

Models that never cross are right-censored at 2021-2023.

## 6. Maintenance policies

Each frozen model is additionally followed under four maintenance policies, each applied at the start of each evaluation cycle using only data available up to and including the previous cycle:

| | Policy |
|---|---|
| P0 | none (frozen) |
| P1 | intercept recalibration |
| P2 | intercept and slope recalibration |
| P3 | recalibration triggered when calibration-in-the-large leaves [−0.2, 0.2] on the log-odds scale |
| P4 | full refit |

Reported: median useful life under each policy, the gain over P0, and the **unrescuable fraction** — the proportion of model–task combinations for which no policy keeps net benefit above the default strategy.

## 7. Analysis

**Primary.** Kaplan–Meier estimate of useful life across the portfolio under P0, with median and interquartile range in years.

**Determinants.** Cox proportional hazards with covariates fixed now: base-period outcome prevalence; absolute change in prevalence per decade across the window; predictor set (Q vs Q+E); model family; base-period effective sample size; base-period AUROC. Proportional hazards checked by scaled Schoenfeld residuals; if violated, an accelerated failure time model with the same covariates is reported instead, and the switch is disclosed.

**Retirement rule comparison.** The net-benefit boundary is compared against three alternatives applied to the same frozen models: AUROC falling more than 0.05 below base period; calibration-in-the-large outside [−0.2, 0.2]; calibration slope outside [0.8, 1.25]. Reported: the cycle each rule fires, the lead or lag in years relative to the net-benefit rule, and the number of cycles each rule would have left a model in service after its net benefit had already fallen to the default strategy.

**Multiplicity.** The primary estimand is the pooled useful-life distribution across the portfolio. Per-task estimates are descriptive and carry no inferential claim; no per-task p-values are reported.

## 8. Survey design and missing data

All estimates are design-based using `WTMEC2YR`, with `SDMVSTRA` as strata and `SDMVPSU` as primary sampling units. Confidence intervals by bootstrap resampling of PSUs within strata, 500 replicates, seed fixed at 20261106.

Missing predictors are imputed once within cycle: median for continuous, mode for binary. Poverty–income ratio, the only predictor missing in more than 5% overall, additionally carries a missing indicator. Observations missing the outcome for a given task are excluded from that task only.

## 9. Measurement-shift handling, fixed now

Blood pressure was measured auscultatorily through 2017-2018 and oscillometrically in 2021-2023. For `y_htn_measured`, the **primary analysis stops at 2017-2018 (ten cycles)**; the 2021-2023 cycle is added in a pre-specified sensitivity analysis and the difference is reported. All other tasks use all eleven cycles.

## 10. What would falsify the paper's claim

If, under P0, useful life is not reached for most of the portfolio within 24 years, the claim that models have a finite and estimable useful life on this substrate is not supported, and we report that. If maintenance (P1–P4) restores net benefit above the default strategy in essentially every case, the "unrescuable fraction" is near zero and the argument for retirement criteria is weakened; we report that too. Neither outcome will be reframed after the fact.
