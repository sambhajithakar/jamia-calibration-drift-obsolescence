#!/usr/bin/env python3
"""
50_prevalence_threshold_signal.py

POST HOC. Written after the stage-1 and stage-2 results were known. This is
hypothesis-generating and is labelled as such wherever it is reported.

Question
--------
The manuscript's Discussion and Future work say that no routine monitoring
signal distinguishes a population that has drifted away from its decision
threshold from a fitted function that has decayed, and that comparing
contemporaneous outcome prevalence with the decision threshold would draw
that distinction. This script tests whether that comparison, computed from
quantities a monitoring system already holds, carries information about
loss of clinical usefulness in this portfolio.

Pre-commitments made BEFORE any operating characteristic was computed
---------------------------------------------------------------------
1. Reference standard: the manuscript's primary (post-hoc) utility rule,
   unchanged, as implemented in 12_boundary_v2.py:
       fail_v2 = nb <= 0.50 * nb_ref
   where nb_ref is the net benefit at the first evaluation cycle under
   policy P0. Combination-level obsolescence requires two consecutive
   breached cycles (K = 2), exactly as in the manuscript. The `fail`
   column of followup_metrics.csv is the PRE-SPECIFIED boundary that the
   manuscript reports as near-vacuous, and is NOT used here.

2. Candidate signal: relative prevalence drift,
       d_c = 1 - prev_c / prev_base
   the proportional fall in design-weighted outcome prevalence since the
   base period. prev_base is the first evaluation cycle, matching how the
   other rules are anchored.

3. We do NOT choose a single cut-off. Choosing one after seeing which value
   separates the known failures would be the same error the manuscript
   criticises in acceptability boundaries. Instead the full operating
   characteristic is reported across all attainable cut-offs, so the reader
   sees whether the signal discriminates at all, and at what cost.

4. Comparators: the three conventional rules exactly as defined in the
   manuscript - calibration-in-the-large outside [-0.20, 0.20], calibration
   slope outside [0.80, 1.25], AUROC at least 0.05 below its first value.

5. Whatever the result, it is reported. If the signal does not discriminate,
   that is the finding.

Outputs
-------
out/prevalence_threshold_signal.csv   per-evaluation-point values
out/prevalence_threshold_roc.csv      operating characteristic across cut-offs
stdout                                summary for the manuscript
"""

import csv
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

CITL_BAND = (-0.20, 0.20)
SLOPE_BAND = (0.80, 1.25)
AUROC_DROP = 0.05


def load_p0():
    rows = []
    with open(os.path.join(OUT, "followup_metrics.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["policy"] != "P0":
                continue
            rows.append(
                dict(
                    task=r["task"],
                    predset=r["predset"],
                    family=r["family"],
                    cycle=r["cycle"],
                    midyear=float(r["midyear"]),
                    prev=float(r["prev"]),
                    thr=float(r["thr"]),
                    auroc=float(r["auroc"]),
                    citl=float(r["citl"]),
                    slope=float(r["slope"]),
                    nb=float(r["nb"]),
                )
            )
    return rows


def first_sustained(series, key, k=2):
    """First midyear at which `key` is True for k consecutive cycles, else None."""
    flags = [bool(r[key]) for r in series]
    for i in range(len(flags) - k + 1):
        if all(flags[i:i + k]):
            return series[i]["midyear"]
    return None


def main():
    rows = load_p0()
    combos = defaultdict(list)
    for r in rows:
        combos[(r["task"], r["predset"], r["family"])].append(r)
    for k in combos:
        combos[k].sort(key=lambda r: r["midyear"])

    recs = []
    for k, series in combos.items():
        base = series[0]
        # manuscript's primary boundary: nb <= 0.50 * nb at first evaluation cycle
        for r in series:
            r["fail"] = r["nb"] <= 0.50 * base["nb"]
        for r in series:
            prev_drift = 1.0 - (r["prev"] / base["prev"]) if base["prev"] > 0 else float("nan")
            gap = (r["thr"] - r["prev"]) / r["thr"] if r["thr"] > 0 else float("nan")
            base_gap = (base["thr"] - base["prev"]) / base["thr"] if base["thr"] > 0 else float("nan")
            recs.append(
                dict(
                    task=k[0], predset=k[1], family=k[2],
                    cycle=r["cycle"], midyear=r["midyear"],
                    prev=r["prev"], prev_base=base["prev"], thr=r["thr"],
                    prev_drift=prev_drift,
                    gap=gap, gap_base=base_gap, gap_widen=gap - base_gap,
                    nb=r["nb"], nb_base=base["nb"], fail=r["fail"],
                    citl_breach=not (CITL_BAND[0] <= r["citl"] <= CITL_BAND[1]),
                    slope_breach=not (SLOPE_BAND[0] <= r["slope"] <= SLOPE_BAND[1]),
                    auroc_breach=r["auroc"] <= base["auroc"] - AUROC_DROP,
                )
            )

    with open(os.path.join(OUT, "prevalence_threshold_signal.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0].keys()))
        w.writeheader()
        w.writerows(recs)

    n = len(recs)
    nfail = sum(r["fail"] for r in recs)
    print(f"evaluation points: {n}   utility failures: {nfail} ({100*nfail/n:.1f}%)")
    print()

    def oc(flag_fn, label):
        tp = sum(1 for r in recs if flag_fn(r) and r["fail"])
        fp = sum(1 for r in recs if flag_fn(r) and not r["fail"])
        fn = sum(1 for r in recs if not flag_fn(r) and r["fail"])
        fired = tp + fp
        ppv = tp / fired if fired else float("nan")
        sens = tp / (tp + fn) if (tp + fn) else float("nan")
        print(f"{label:<34} fired {fired:>4}  TP {tp:>3}  FP {fp:>4}  PPV {ppv:5.3f}  sens {sens:5.3f}")
        return dict(rule=label, fired=fired, tp=tp, fp=fp, ppv=ppv, sens=sens)

    print("CONVENTIONAL RULES (as in the manuscript)")
    comp = [
        oc(lambda r: r["citl_breach"], "calibration-in-the-large"),
        oc(lambda r: r["slope_breach"], "calibration slope"),
        oc(lambda r: r["auroc_breach"], "AUROC drop >= 0.05"),
    ]
    print()
    print("CANDIDATE SIGNAL: relative prevalence drift, full operating characteristic")
    print("(cut-off d: flag when prevalence has fallen by at least d since base period)")
    print()
    cuts = [i / 100 for i in range(0, 51, 2)]
    roc = []
    for d in cuts:
        tp = sum(1 for r in recs if r["prev_drift"] >= d and r["fail"])
        fp = sum(1 for r in recs if r["prev_drift"] >= d and not r["fail"])
        fn = sum(1 for r in recs if not (r["prev_drift"] >= d) and r["fail"])
        fired = tp + fp
        ppv = tp / fired if fired else float("nan")
        sens = tp / (tp + fn) if (tp + fn) else float("nan")
        roc.append(dict(cutoff=d, fired=fired, tp=tp, fp=fp, ppv=ppv, sens=sens))
        print(f"  d >= {d:4.2f}   fired {fired:>4}  TP {tp:>3}  FP {fp:>4}  PPV {ppv:5.3f}  sens {sens:5.3f}")

    with open(os.path.join(OUT, "prevalence_threshold_roc.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["cutoff", "fired", "tp", "fp", "ppv", "sens"])
        w.writeheader()
        w.writerows(roc)

    # ---- combination level, two consecutive breached cycles, as in the manuscript ----
    print()
    print("COMBINATION LEVEL (K=2 consecutive cycles, as in the manuscript)")
    by_combo = defaultdict(list)
    for r in recs:
        by_combo[(r["task"], r["predset"], r["family"])].append(r)
    for kk in by_combo:
        by_combo[kk].sort(key=lambda x: x["midyear"])

    def combo_oc(flagkey, label, extra=None):
        tp = fp = fn = 0
        for kk, s in by_combo.items():
            if extra is not None:
                for r in s:
                    r["_tmp"] = extra(r)
                fired = first_sustained(s, "_tmp") is not None
            else:
                fired = first_sustained(s, flagkey) is not None
            failed = first_sustained(s, "fail") is not None
            if fired and failed: tp += 1
            elif fired and not failed: fp += 1
            elif failed and not fired: fn += 1
        n_fired = tp + fp
        ppv = tp / n_fired if n_fired else float("nan")
        sens = tp / (tp + fn) if (tp + fn) else float("nan")
        print(f"  {label:<40} fired {n_fired:>2}/36  TP {tp:>2}  FP {fp:>2}  PPV {ppv:5.3f}  sens {sens:5.3f}")

    n_obsolete = sum(1 for kk, s in by_combo.items() if first_sustained(s, "fail") is not None)
    print(f"  combinations obsolete by the utility rule: {n_obsolete}/36")
    combo_oc("citl_breach", "calibration-in-the-large")
    combo_oc("slope_breach", "calibration slope")
    combo_oc("auroc_breach", "AUROC drop >= 0.05")
    for d in (0.05, 0.10, 0.15, 0.20):
        combo_oc(None, f"prevalence drift >= {d:.2f}", extra=lambda r, d=d: r["prev_drift"] >= d)

    # ---- the manuscript's actual claim: does prevalence drift SEPARATE the
    #      calibration-in-the-large true alarms from its false alarms? ----
    print()
    print("THE MANUSCRIPT'S CLAIM, TESTED DIRECTLY")
    print("Among combinations where calibration-in-the-large fired, does prevalence")
    print("drift distinguish the genuine utility failures from the false alarms?")
    print()
    hdr = f"  {'combination':<34} {'CITL':>5} {'failed':>7} {'max drift':>10} {'min prev/base':>14}"
    print(hdr)
    rowsout = []
    for kk, s in by_combo.items():
        citl_fired = first_sustained(s, "citl_breach") is not None
        failed = first_sustained(s, "fail") is not None
        if not citl_fired:
            continue
        max_drift = max(r["prev_drift"] for r in s)
        min_ratio = min(r["prev"] / r["prev_base"] for r in s)
        name = f"{kk[0].replace('y_','')}/{kk[1]}/{kk[2]}"
        rowsout.append((name, failed, max_drift, min_ratio))
    for name, failed, md, mr in sorted(rowsout, key=lambda x: -x[2]):
        print(f"  {name:<34} {'yes':>5} {('YES' if failed else 'no'):>7} {md:>10.3f} {mr:>14.3f}")

    tf = [r for r in rowsout if r[1]]
    ff = [r for r in rowsout if not r[1]]
    if tf and ff:
        print()
        print(f"  true alarms  (n={len(tf)}): max drift {min(r[2] for r in tf):+.3f} to {max(r[2] for r in tf):+.3f}")
        print(f"  false alarms (n={len(ff)}): max drift {min(r[2] for r in ff):+.3f} to {max(r[2] for r in ff):+.3f}")
        sep = min(r[2] for r in tf) > max(r[2] for r in ff)
        print(f"  cleanly separable by prevalence drift alone: {sep}")

    # ---- second pre-specified candidate: the threshold-prevalence GAP itself ----
    print()
    print("SECOND PRE-SPECIFIED CANDIDATE: the threshold-prevalence gap")
    print("gap = (threshold - contemporaneous prevalence) / threshold")
    print()
    print(f"  {'combination':<34} {'failed':>7} {'gap base':>9} {'gap max':>9} {'widening':>9}")
    g_rows = []
    for kk, s in by_combo.items():
        failed = first_sustained(s, "fail") is not None
        gmax = max(r["gap"] for r in s)
        gbase = s[0]["gap_base"]
        g_rows.append((f"{kk[0].replace('y_','')}/{kk[1]}/{kk[2]}", failed, gbase, gmax, gmax - gbase))
    for name, failed, gb, gm, gw in sorted(g_rows, key=lambda x: -x[3]):
        mark = "YES" if failed else "no"
        print(f"  {name:<34} {mark:>7} {gb:>9.3f} {gm:>9.3f} {gw:>9.3f}")

    tf = [r for r in g_rows if r[1]]
    ff = [r for r in g_rows if not r[1]]
    print()
    for idx, lab in ((3, "gap (max)"), (4, "gap widening")):
        lo_t, hi_t = min(r[idx] for r in tf), max(r[idx] for r in tf)
        lo_f, hi_f = min(r[idx] for r in ff), max(r[idx] for r in ff)
        sep = lo_t > hi_f
        print(f"  {lab:<14} true {lo_t:+.3f}..{hi_t:+.3f}   false {lo_f:+.3f}..{hi_f:+.3f}   separable: {sep}")

    print()
    print("Which combinations ever failed, and their prevalence drift at failure:")
    seen = set()
    for r in sorted(recs, key=lambda x: (x["task"], x["predset"], x["family"], x["midyear"])):
        if r["fail"]:
            k = (r["task"], r["predset"], r["family"])
            if k in seen:
                continue
            seen.add(k)
            print(f"  {r['task']:<20} {r['predset']}/{r['family']:<4} {r['cycle']}  "
                  f"prev {r['prev_base']:.3f} -> {r['prev']:.3f}  drift {r['prev_drift']:+.3f}  "
                  f"thr {r['thr']:.2f}  NB {r['nb_base']:+.4f} -> {r['nb']:+.4f}")


if __name__ == "__main__":
    main()
