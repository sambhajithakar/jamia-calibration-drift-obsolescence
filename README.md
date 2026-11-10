# Calibration drift is not obsolescence — analysis archive

Code, specifications and derived data for:

> Thakar SB, et al. *Calibration drift is not obsolescence: a 20-year decision-analytic
> study of 36 frozen clinical prediction models.* Submitted to the Journal of the
> American Medical Informatics Association, November 2026.

Everything reported in the manuscript is reproducible from this archive.

## Source data

Raw data are the public NHANES Continuous files for the eleven cycles 1999-2000 through
2021-2023 (the 2019-2020 cycle was not released as a standalone public file). They are
**not redistributed here**; download them from
https://wwwn.cdc.gov/nchs/nhanes/ . `file_inventory.csv` lists all 77 `.xpt` files used,
with the SHA-256 of each, so an independent analyst can confirm they have the same files.

## Order of execution

| Script | What it does | Writes |
|---|---|---|
| `01_harmonisation_audit.py` | inventories the 77 `.xpt` files, checks variable coverage across cycles | `file_inventory.csv`, `variable_coverage.csv` |
| `02_build_analytic.py` | recode and merge; **computes no predictor–outcome association** | `analytic.csv.gz` |
| `10_fit_and_follow.py` | fits the frozen models, follows them across evaluation cycles under five maintenance policies | `followup_metrics.csv`, `task_eligibility.csv` |
| `11_survival.py` | executes PRESPECIFICATION_v1 §5–7 (the pre-specified analysis that falsified its own hypothesis) | `V1_PRESPECIFIED_RESULTS.txt` |
| `12_boundary_v2.py` | re-derives event times under the post-hoc boundary, without refitting | console output |
| `31_threshold_mechanism.py` | executes PRESPECIFICATION_v1 §4 (thresholds at 0.5x and 2x) and §9 (hypertension with 2021-2023 included), and tests the mechanism directly by holding the frozen predictions fixed while the threshold tracks contemporaneous prevalence | `out/threshold_mechanism.csv` |
| `40_bootstrap_ci.py` | executes PRESPECIFICATION_v1 §8: design-based 95% confidence intervals by resampling PSUs within strata, 500 replicates, seed 20261106. Asserts that its point estimates reproduce `followup_metrics.csv` exactly | `out/bootstrap_ci.csv` |
| `20_tables.py` | manuscript Tables 1–4 | `tables/*.csv`, `tables/key_numbers.json` |
| `30_fit_weighted_threshold.py` | sensitivity analysis: identical to `10_fit_and_follow.py` except that the decision threshold is the **design-based** base-period prevalence rather than the unweighted sample proportion | `out_wtd/followup_metrics.csv`, `out_wtd/task_eligibility.csv` |
| `21_figures.py` | manuscript Figures 1–6 | `figures/*.png` |

Scripts are numbered so that the order in which they were written is visible: everything
below `10_` predates the pre-specification and is structural only.

## Specifications

- `PRESPECIFICATION_v1.md` — written before any predictor–outcome association was computed.
  SHA-256 `5daf3197b26be2714c3f5f2b6fc240873792ce5f487360869f677e598af4aac8`.
  Its §10 falsification conditions **were met**; see `RESULTS_MEMO_01.md`.
- `PRESPECIFICATION_v2_posthoc.md` — the replacement boundary, written after seeing the v1
  result and labelled post-hoc throughout.
  SHA-256 `999a05d1a263f59d3150e40fa9cc3c49258b0f16e1d9cbeecdb39d522a9e1777`.
- `RESULTS_MEMO_01.md` — the contemporaneous record of the v1 failure, its diagnosis, and
  the options considered before any re-analysis was run.

Verify with `sha256sum -c SHA256SUMS`.

## Corrections made during preparation

Two errors were found and fixed before any result was reported, and both are recorded here
rather than silently corrected:

1. **Intercept-only recalibration (policy P1) was implemented incorrectly.** The intercept
   from a two-parameter recalibration was paired with a slope fixed at 1.0, which is not a
   valid intercept-only fit. Replaced with a weighted offset fit
   (`calib_intercept_only` in `10_fit_and_follow.py`). Before the fix P1 appeared to make
   models worse (14/36); after it, P1 gives 4/36, identical to P0, P2 and P3. The erroneous
   result was never reported.
2. **The secondary utility boundary was coded differently from the specification.** An early
   draft of `20_tables.py` used "net benefit falls by ≥1 net true positive per 1000 from its
   own initial value"; `PRESPECIFICATION_v2_posthoc.md` specifies "net benefit adds <1 net
   true positive per 1000 over the better default strategy". The script now implements the
   specification.
3. **The survey strata and primary sampling units were read and never used.** The
   Methods described every quantity as design-based, but `SDMVSTRA` and `SDMVPSU` were
   assigned to local variables and never referenced again. Closed by running the
   specified bootstrap rather than by softening the claim — see the deviations section
   below.
4. **A claim in the manuscript was false and has been withdrawn.** An earlier draft said
   that the threshold–prevalence failure was invisible to any statistic computed from a
   model's own predictions. Calibration-in-the-large is such a statistic, and by these
   data it fired on the affected models. The claim is replaced by the narrower and true
   one: a band breach cannot distinguish a drifted population from a decayed model, and
   the two call for opposite responses.

## Sensitivity analysis: how the decision threshold was derived

Thresholds in the primary analysis are the **unweighted** base-period sample proportion,
rounded to 0.01, while every evaluation quantity (calibration, net benefit, prevalence) is
**design-based** (`WTMEC2YR`, `SDMVSTRA`, `SDMVPSU`). Because NHANES oversamples, the two
estimators differ: the threshold sits above the design-based base-period prevalence for all
ten tasks, by 0.007 to 0.070. `30_fit_weighted_threshold.py` repeats the whole follow-up
with design-based thresholds, changing nothing else.

Result: every reported quantity is unchanged except the utility breach rate, which falls
from 6.6% to 4.1% of evaluation points. Trigger counts (23/20/4/0), false alarms (19 for
calibration-in-the-large, 17 for slope), positive predictive values (0.17 and 0.15), the
single failing task, and the 100%/89%/44%/36% apparent survival are identical. The
signal–utility discordance therefore does not depend on how the threshold was derived.

## Deviations from the pre-specification, and how they were closed

Three items in `PRESPECIFICATION_v1.md` were not executed in the first analysis pass and
were found by an independent audit of this archive against the manuscript. All three have
since been run; none changed a reported conclusion.

1. **§8, design-based confidence intervals.** The first pass read `SDMVSTRA` and `SDMVPSU`
   and then used neither: point estimates were weight-only and no interval was computed.
   `40_bootstrap_ci.py` now performs the specified PSU-within-stratum bootstrap (500
   replicates, seed 20261106) and asserts that its point estimates reproduce the primary
   run to within 1e-15. Result: of the breaches, the whole 95% interval lies outside the
   band at 33.4% of evaluation points for calibration slope and 25.9% for
   calibration-in-the-large, against raw breach rates of 49.4% and 45.3% — about a third of
   calibration alarms are not separable from sampling error.
2. **§4, thresholds at 0.5x and 2x.** Now reported: 0 of 36 combinations obsolete at half
   the threshold, 4 of 36 at the primary value, 11 of 36 at double.
3. **§9, hypertension including 2021-2023.** Now reported: all four variants obsolete either
   way, so the pre-specified stop at 2017-2018 does not drive the result.

`31_threshold_mechanism.py` additionally runs a test that is not in either specification and
is labelled post-hoc in the manuscript: holding the frozen predictions fixed and letting only
the decision threshold track contemporaneous weighted prevalence. Obsolescence falls from 4 of
36 to 1 of 36, and the hypertension models retain 72-82% of their initial net benefit at the
2009-2010 trough instead of -17% to +11%. This is the direct evidence for the manuscript's
claim that the loss was in the threshold rather than the model.

## The mechanism test

`31_threshold_mechanism.py` runs a test that is in neither specification and is labelled
post-hoc in the manuscript: the frozen predictions are held exactly as they are and only the
decision threshold is allowed to track contemporaneous weighted prevalence. **All four
hypertension failures disappear** — one obesity model crosses the boundary that had not
before, so the count falls from 4/36 to 1/36 with a different member — and at the 2009–2010
trough the hypertension models retain 72–82% of their initial net benefit instead of −17% to
+11%. Nothing about the model changed; only the threshold did. This is the direct evidence
for the manuscript's claim that the loss was in the threshold rather than the model.

## Environment

The exact versions the reported results were produced with:

| | |
|---|---|
| Python | 3.13.18 |
| numpy | 2.5.3 |
| scipy | 1.18.1 |
| scikit-learn | 1.9.1 |
| pandas | 2.3.3 |
| lifelines | 0.30.3 |
| matplotlib | 3.11.2 |
| pyreadstat | 1.3.7 |

scikit-learn's estimators are version-sensitive: a different minor version may give
fitted coefficients that differ in the last decimal places, and therefore metrics that
differ in the last decimal places. Nothing in the manuscript turns on a difference of
that size, but a reproducer who wants bit-identical output should pin these versions.

## Licence

Code: MIT. Derived data: CC0, consistent with the public domain status of the NHANES source.

---

## Version history

**1.2.0 (11 November 2026).** Adds the post-hoc test of prevalence drift and the
threshold-prevalence gap as a candidate monitoring signal (`50_prevalence_threshold_signal.py`,
`out/prevalence_threshold_signal.csv`, `out/prevalence_threshold_roc.csv`,
`tables/table5_threshold_signal.csv`). The script states its definitions and its reporting
rule in its docstring before any operating characteristic is computed, and the result is
negative: no cut-off on either quantity separates the four combinations that lost half their
net benefit from the 32 that did not. The reference standard is the manuscript's primary
(post-hoc) utility rule, nb <= 0.50 * nb_ref with two consecutive breached cycles, as
implemented in 12_boundary_v2.py; the near-vacuous pre-specified boundary is not used here.
No existing script, output or point estimate is changed by this version.

**1.1.1 (8 November 2026).** Three code-level corrections found by a pre-submission audit
of this archive against the manuscript. None changes a point estimate.

1. `40_bootstrap_ci.py` derived its per-cell seed from `hash()` of a tuple of strings.
   Python salts `hash()` per process unless `PYTHONHASHSEED` is set, so the replicate draws
   were **not** reproducible across runs — contrary to what PRESPECIFICATION_v1 §8 and the
   manuscript both claim. The seed is now a SHA-256 digest of the cell key, and the whole
   bootstrap has been re-run. The end-of-run assertion against `followup_metrics.csv`
   confirms the point estimates are unchanged. Two reported percentages moved:

   | | 1.1.0 | 1.1.1 |
   |---|---|---|
   | whole 95% interval outside the band, calibration slope | 34.4% | **33.4%** |
   | whole 95% interval outside the band, calibration-in-the-large | 25.6% | **25.9%** |
   | slope breaches not separable from sampling error | 30% | **32%** |
   | calibration-in-the-large breaches not separable | 43% | 43% |

2. `20_tables.py` bound `tot` to the 320 evaluation points and then rebound it to the 36
   model–tasks inside a later loop, so `tables/key_numbers.json` recorded
   `n_eval_points: 36`. Fixed at the source rather than in the JSON. The typographic
   normalisation of the delivered tables moved into the same script, so the tables in this
   archive are now produced by it rather than edited afterwards.

3. `21_figures.py` — Figure 1 stated that the four metrics were "design-based (NHANES MEC
   weights, strata, PSUs)". They are not: the point estimates use `WTMEC2YR` only, and the
   strata and PSUs enter the bootstrap alone. This had been corrected in the manuscript's
   Methods but left standing in the figure. Regenerating the figures produced the other
   five byte-identical.

Also added: `PRESPECIFICATION_v1.sha256`, which was missing; `STATISTICAL_INTEGRITY.md`,
which records every error found by verification and its fix; and exact library versions
above. `SHA256SUMS` regenerated over all 46 files.

**1.1.0 (7 November 2026).** Closed three deviations from the pre-specification: the
bootstrap confidence intervals of §8, and the threshold sensitivity analyses of §4 and §9.
Added the post-hoc mechanism test and the weighted-threshold sensitivity analysis.

**1.0.0 (6 November 2026).** First deposit.
