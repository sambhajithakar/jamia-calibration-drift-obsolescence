# Pre-specification v2 — POST-HOC, written after seeing the v1 result
6 November 2026.

## Status, stated plainly
This is **not** a blind pre-specification. It was written after the v1 analysis
(sha256 5daf3197…) returned its result. Every number produced under v2 must be
reported as post-hoc. v1 remains the primary analysis and its result stands as
reported in RESULTS_MEMO_01.md.

## Why v2 exists
The v1 failure boundary set the decision threshold to the base-period prevalence.
At that threshold the net benefit of treat-all is exactly zero by construction, and
because prevalence fell over the window it became negative at 66% of evaluation
points. The boundary therefore collapsed to "the model must be worse than doing
nothing", a test of active harm rather than of lost usefulness. It is too lenient
to answer the question the study asks.

## The v2 boundary
For each model-task combination let NB_ref be the model's net benefit at its **first
evaluation cycle** (2003-2004) — out of sample, never the training data.

A model is **obsolete** at the first cycle at which

        NB_model(cycle)  <=  0.50 x NB_ref

sustained for k = 2 consecutive evaluation cycles. In words: the model has lost half
the clinical value it was delivering when it entered service.

Sensitivity: retention fractions 0.25 and 0.75; k = 1 and 3.

Secondary absolute boundary: NB_model falls below NB_default + 0.001, i.e. the model
adds fewer than 1 net true positive per 1000 people screened.

## What is unchanged from v1
Portfolio, base period, models, maintenance policies, survey design handling,
imputation, the measurement-shift rule for blood pressure, and the comparison
against AUROC- and calibration-based rules. No model is refitted; v2 re-derives
event times from the metrics already computed under v1.

## What would still falsify the argument
If under this boundary the calibration rules no longer fire far more often than the
utility rule, then calibration is not a poor proxy for clinical usefulness, the
surviving finding from v1 does not hold, and we report that and stop.
