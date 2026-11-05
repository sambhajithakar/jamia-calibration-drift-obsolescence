# NHANES harmonisation audit — result
Run 5 November 2026 on a frozen, SHA-256-checksummed copy of 77 `.xpt` files.
Script: `03_analysis/01_harmonisation_audit.py`. Outputs: `out/file_inventory.csv`, `out/variable_coverage.csv`.
This audit inspected file structure and variable names only. No outcome relationship was examined, so the pre-specification that follows is still blind.

## Headline

**The study window is the full 11 cycles, 1999–2000 through 2021–2023 — about 24 years of calendar time.** That is roughly double what v2 of the protocol assumed, and more than twice the span of the nearest comparable work (Davis et al., JAMIA 2025, used 11 years of VA data).

## Cycles recovered

1999-2000, 2001-2002, 2003-2004, 2005-2006, 2007-2008, 2009-2010, 2011-2012, 2013-2014, 2015-2016, 2017-2018, 2021-2023.

**2019–2020 is absent by design**, not by omission: NHANES could not complete that cycle because of COVID-19 and released it only combined with 2017–2018 as the "pre-pandemic" P_ files. The gap must be stated in the paper and handled as a gap in the time axis, not silently closed. It is also substantively interesting: the obsolescence clock keeps running through a period with no observation, which is exactly what happens to a deployed model during a service interruption.

## Row counts per cycle (range across cycles)

| Component | Cycles | Rows (min–max) |
|---|---|---|
| DEMO (demographics) | 11 | 9,254 – 11,933 |
| BMX (body measures) | 11 | 8,704 – 10,477 |
| DIQ (diabetes questionnaire) | 11 | 8,897 – 11,744 |
| MCQ (medical conditions) | 11 | 8,897 – 11,744 |
| SMQ (smoking) | 11 | 4,880 – 9,015 |
| BPX (blood pressure, auscultatory) | 10 | 8,704 – 10,477 |
| BPXO (blood pressure, oscillometric) | 1 (2021-2023) | 7,801 |
| GHB (HbA1c) | 8 (2005 onward) | 6,401 – 7,369 |
| LAB10 / L10 (HbA1c, earlier name) | 3 (1999–2004) | 6,758 – 7,445 |

## Variables available in **every** cycle

- **DEMO (17)** including the full survey design set — `SDMVPSU`, `SDMVSTRA`, `WTMEC2YR`, `WTINT2YR` — plus `RIDAGEYR`, `RIAGENDR`, `RIDRETH1`, `DMDEDUC2`, `INDFMPIR`, `RIDEXPRG`, `SDDSRVYR`.
  Design-based analysis is therefore possible across the whole window, which matters: a survey-naive analysis would not survive JAMIA review.
- **BMX (18)** including `BMXBMI`, `BMXWAIST`, `BMXHT`, `BMXWT`, `BMXARMC`.
- **MCQ (18)** including `MCQ160A–F` (arthritis, heart failure, coronary heart disease, angina, heart attack, stroke) and `MCQ220` (cancer).
- **SMQ (4)** including `SMQ020` (ever smoked 100 cigarettes) and `SMQ040` (currently smokes).
- **DIQ (4)** including `DIQ010` (ever told by a doctor you have diabetes).

## The two measurement breaks, and what they mean

**1. HbA1c is available in all 11 cycles — but under three different file names.**
`LBXGH` appears in `LAB10` (1999–2000), `L10` (2001–2004) and `GHB` (2005–2023). The variable name is stable; only the container changed. So the **primary outcome can be defined across the full 24 years**.

**2. Blood pressure changed measurement method in 2021–2023.**
Cycles 1999–2018 use auscultatory readings (`BPXSY1–4`, `BPXDI1–4`). The 2021–2023 cycle uses oscillometric readings (`BPXOSY1–3`, `BPXODI1–3`) — a different instrument, different variable names, and a known systematic difference between methods.

This is not a nuisance. It is the cleanest available example of the thing the paper is about: an environmental change that a deployed model's monitoring would never see, because nothing about the model changed. Two options, and the choice must be pre-specified:
- (a) restrict blood-pressure-dependent tasks to 1999–2018 (10 cycles) and report the HbA1c-based tasks over the full 11; or
- (b) carry all tasks to 11 cycles and treat 2021–2023 as a declared measurement-shift cycle, analysed both with and without.

**Recommendation: (a) as primary, (b) as a pre-specified sensitivity analysis**, with the discrepancy reported either way. Hiding it would be the error; using it is the contribution.

## What this means for the design

- The primary prediction task — **undiagnosed diabetes**, defined as HbA1c ≥ 6.5% in an adult who answers no to `DIQ010` — is feasible over all 11 cycles. Prevalence, screening practice, obesity and treatment all changed substantially over this period, so the drift is real and documented rather than simulated.
- A portfolio of secondary tasks is feasible from `MCQ160A–F`, `MCQ220` and anthropometric outcomes, all with complete cycle coverage.
- Predictor classes separate cleanly, which the design needs as an analysable covariate: questionnaire only (DEMO + SMQ + MCQ), plus examination (BMX), plus laboratory (HbA1c, lipids).
- Survey design variables are complete, so every estimate can be design-based with bootstrap confidence intervals over PSUs within strata.

## Risk retired

The largest risk named in protocol v2 — that harmonisation across cycles would force the window to be trimmed — **did not materialise**. The window is longer than planned, not shorter. The remaining harmonisation work is recoding, not reconciliation.
