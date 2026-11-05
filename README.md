# Calibration drift is not obsolescence

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23988407.svg)](https://doi.org/10.5281/zenodo.23988407)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Analysis code, pre-specifications and derived data for:

> Thakar SB, Ghorpade PN. *Calibration drift is not obsolescence: a 20-year
> decision-analytic study of 36 frozen clinical prediction models.* Submitted to the
> Journal of the American Medical Informatics Association, November 2026.

**Archived release with DOI:** [10.5281/zenodo.23988407](https://doi.org/10.5281/zenodo.23988407)

## What the study found

Thirty-six model-task combinations were trained on the 1999-2002 cycles of NHANES,
frozen, and followed across nine later cycles to 2021-2023 - 20.5 years of calendar
time, 320 model-cycle evaluation points.

| Monitoring signal | Evaluation points breached | Model-tasks triggered |
|---|---|---|
| Calibration slope outside [0.80, 1.25] | 158/320 (49.4%) | 20/36 |
| Calibration-in-the-large outside [-0.20, 0.20] | 145/320 (45.3%) | 23/36 |
| Net benefit fell to half its initial value | 21/320 (6.6%) | 4/36 |
| AUROC fell >0.05 below its initial value | 0/320 (0.0%) | 0/36 |

Among the 32 combinations that never lost half their clinical value,
calibration-in-the-large fired on 19 and calibration slope on 17. Applied to the same
portfolio, the four rules give four different answers to "how many of these models are
still fit for use at 20 years?" - 100%, 89%, 44% and 36%.

Everything reported in the manuscript is reproducible from this repository.

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
| `20_tables.py` | manuscript Tables 1–4 | `tables/*.csv`, `tables/key_numbers.json` |
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

## Environment

Python 3.13; pandas, numpy, scipy, scikit-learn, lifelines, pyreadstat, matplotlib.

## Licence

Code: MIT. Derived data: CC0, consistent with the public domain status of the NHANES source.
