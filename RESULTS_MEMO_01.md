# Results memo 1 — the pre-specified analysis contradicted the hypothesis
6 November 2026. Reports the output of `10_fit_and_follow.py` and `11_survival.py` executing PRESPECIFICATION_v1.md (sha256 5daf3197…).

## What was expected

That frozen models would lose clinical usefulness over 24 years; that maintenance would postpone but not always prevent this; and that calibration would look acceptable while net benefit was already gone — so that calibration-based monitoring would miss the moment a model should be retired.

## What actually happened

**Portfolio:** 36 model–task combinations (10 tasks, Q and Q+E predictor sets, logistic regression and gradient-boosted trees), frozen on 1999–2002 and followed to 2021–2023.

| | |
|---|---|
| Obsolete within the window (P0, frozen) | **3 / 36 (8.3%)** |
| Median useful life | **not reached** |
| Survival at 20 years | 0.92 |
| Unrescuable under all five maintenance policies | **0 / 36 (0%)** |

Maintenance made it rarer still: P1 2/36, P2 1/36, P3 1/36, P4 (full refit) 0/36.

**Both falsification conditions written into §10 of the pre-specification were met.** The claim that these models have a finite, estimable useful life on this substrate is not supported, and the "unrescuable fraction" argument for retirement criteria is not supported. Per the pre-specification, this is reported, not reframed.

## Why — and this part is a flaw in my own pre-specification

The pre-specified failure boundary was: net benefit ≤ max(NB of treat-all, 0), with the decision threshold *t* set to the base-period prevalence.

At *t* = prevalence, the net benefit of treat-all is **exactly zero by construction**:

    NB_all = prev − (1 − prev) × t/(1 − t),  and at t = prev this is prev − prev = 0.

Because prevalence fell over the window for most tasks, NB_all went *negative* — it is ≤ 0 at 66% of all evaluation points. So `max(NB_all, 0)` collapsed to 0, and the boundary became "**the model must be worse than doing nothing at all**". That is a test of active harm, not of lost usefulness. It is a very low bar, and models cleared it at 95% of evaluation points.

I chose that threshold blind, in good faith, as the conventional default. It was the wrong choice, and the primary analysis is correspondingly near-vacuous. That has to be said plainly rather than patched over.

## The finding that does survive, and is the more interesting one

Per evaluation point, across the whole portfolio:

| Monitoring signal | Breach rate |
|---|---|
| Calibration-in-the-large outside [−0.2, 0.2] | **45.3%** |
| Calibration slope outside [0.8, 1.25] | **49.4%** |
| Net benefit at or below the default strategy | **4.7%** |

And for the three models that did lose net benefit, the calibration rules fired **4 to 6 years earlier**, not later.

So the direction is the opposite of the hypothesis. Calibration drifted constantly while clinical usefulness almost never went away. **Calibration-based monitoring, used as a trigger to update or retire, would have fired on roughly half of all model-years in which the model was still clinically useful.**

That is a real result, it is policy-relevant, and it speaks directly to an open question Davis, Walsh and Matheny (Front Digit Health 2022) posed and nobody has answered: *"Guidance on whether to update models when statistically significant improvement is possible but updating would not provide a clinically meaningful improvement."*

It also cuts against the working assumption of the whole algorithmovigilance programme — that drift detection is the right trigger for maintenance. On this evidence, drift detection mostly detects drift that does not matter.

## Three honest options

**A. Report v1 as primary, pivot the paper to the surviving finding.**
Keep the pre-specified result exactly as it came out, including the admission that the boundary was too lenient and why. Make the paper about calibration being a poor proxy for clinical usefulness, with the retirement question as the motivation rather than the answer. The pre-specification, including its own falsification clause, becomes a strength: it is a paper about monitoring criteria that reports its own criterion failing.

**B. Write pre-specification v2 with a defensible boundary and re-run, labelling everything from v2 onward as post-hoc.**
A better boundary would be a *margin* — for example net benefit failing to exceed the better default by some minimum clinically worthwhile amount, or expressed as net true positives per 1,000 screened. Both v1 and v2 results get reported, with v2 clearly marked as not pre-specified.

**C. Both: A as the paper, B as a pre-specified sensitivity analysis in the supplement.**

## Recommendation

**C.** The surviving finding is stronger than the one I set out to show, and it is genuinely novel. But the boundary problem must be visible in the main text, not buried — partly because it is honest, and partly because a paper arguing that monitoring criteria are chosen carelessly is in no position to hide a carelessly chosen criterion.

## What I have not done

I have not re-run anything with a different boundary, and I will not until you choose. No result in this memo has been selected after seeing others; everything here is the complete output of the two pre-specified scripts.
