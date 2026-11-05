"""
Build the harmonised analytic dataset across 11 NHANES cycles (1999-2000 .. 2021-2023).

Scope of this script: recoding and merging only. It writes one row per examined adult
per cycle with harmonised predictors, candidate outcomes and survey design variables.
It deliberately computes NO association between any predictor and any outcome, so the
pre-specification written after it is still blind to the results.
"""
import glob, os, hashlib
import numpy as np
import pandas as pd

SRC = "/mnt/user-data/uploads/Downloads/NHANES_JAMIA/_frozen_v1"
OUT = "/home/claude/jamia2026/03_analysis/out"
os.makedirs(OUT, exist_ok=True)

CYCLE = {"": "1999-2000", "_B": "2001-2002", "_C": "2003-2004", "_D": "2005-2006",
         "_E": "2007-2008", "_F": "2009-2010", "_G": "2011-2012", "_H": "2013-2014",
         "_I": "2015-2016", "_J": "2017-2018", "_L": "2021-2023"}
# midpoint year, used as the calendar axis
MIDYEAR = {"1999-2000": 1999.5, "2001-2002": 2001.5, "2003-2004": 2003.5,
           "2005-2006": 2005.5, "2007-2008": 2007.5, "2009-2010": 2009.5,
           "2011-2012": 2011.5, "2013-2014": 2013.5, "2015-2016": 2015.5,
           "2017-2018": 2017.5, "2021-2023": 2022.0}
SUFFIXES = ["", "_B", "_C", "_D", "_E", "_F", "_G", "_H", "_I", "_J", "_L"]

# HbA1c lives in a differently named file in the first three cycles
A1C_FILE = {"": "LAB10", "_B": "L10_B", "_C": "L10_C"}


def read(stem):
    p = os.path.join(SRC, f"{stem}.xpt")
    if not os.path.exists(p):
        return None
    return pd.read_sas(p, format="xport", encoding="latin1")


def keep(df, cols):
    if df is None:
        return None
    have = [c for c in cols if c in df.columns]
    return df[["SEQN"] + [c for c in have if c != "SEQN"]].copy()


frames = []
for suf in SUFFIXES:
    cyc = CYCLE[suf]

    demo = keep(read(f"DEMO{suf}"), ["SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1",
                                     "DMDEDUC2", "INDFMPIR", "RIDEXPRG", "RIDSTATR",
                                     "SDMVPSU", "SDMVSTRA", "WTMEC2YR"])
    bmx  = keep(read(f"BMX{suf}"),  ["SEQN", "BMXBMI", "BMXWAIST", "BMXHT", "BMXWT", "BMXARMC"])
    mcq  = keep(read(f"MCQ{suf}"),  ["SEQN", "MCQ160A", "MCQ160B", "MCQ160C",
                                     "MCQ160E", "MCQ160F", "MCQ220"])
    smq  = keep(read(f"SMQ{suf}"),  ["SEQN", "SMQ020", "SMQ040"])
    diq  = keep(read(f"DIQ{suf}"),  ["SEQN", "DIQ010"])
    a1c  = keep(read(A1C_FILE.get(suf, f"GHB{suf}")), ["SEQN", "LBXGH"])

    # blood pressure: auscultatory to 2017-2018, oscillometric in 2021-2023
    if suf == "_L":
        bp = read("BPXO_L")
        if bp is not None:
            sy = [c for c in ["BPXOSY1", "BPXOSY2", "BPXOSY3"] if c in bp.columns]
            di = [c for c in ["BPXODI1", "BPXODI2", "BPXODI3"] if c in bp.columns]
            bp = pd.DataFrame({"SEQN": bp.SEQN,
                               "SBP": bp[sy].replace(0, np.nan).mean(axis=1),
                               "DBP": bp[di].replace(0, np.nan).mean(axis=1)})
        bp_method = "oscillometric"
    else:
        bpr = read(f"BPX{suf}")
        bp = None
        if bpr is not None:
            sy = [c for c in ["BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4"] if c in bpr.columns]
            di = [c for c in ["BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4"] if c in bpr.columns]
            bp = pd.DataFrame({"SEQN": bpr.SEQN,
                               "SBP": bpr[sy].replace(0, np.nan).mean(axis=1),
                               "DBP": bpr[di].replace(0, np.nan).mean(axis=1)})
        bp_method = "auscultatory"

    d = demo
    for part in [bmx, mcq, smq, diq, a1c, bp]:
        if part is not None:
            d = d.merge(part, on="SEQN", how="left")

    d["cycle"] = cyc
    d["midyear"] = MIDYEAR[cyc]
    d["bp_method"] = bp_method
    frames.append(d)

df = pd.concat(frames, ignore_index=True, sort=False)

# ---------- harmonised recoding ----------
yn = lambda s: s.map({1: 1, 2: 0})          # 1=yes 2=no, everything else missing

df["age"]      = df.RIDAGEYR
df["female"]   = (df.RIAGENDR == 2).astype(float)
df["race"]     = df.RIDRETH1
df["educ_hs"]  = df.DMDEDUC2.where(df.DMDEDUC2 <= 5).ge(3).astype(float)  # >= high school
df["pir"]      = df.INDFMPIR
df["bmi"]      = df.BMXBMI
df["waist"]    = df.BMXWAIST
df["height"]   = df.BMXHT
df["weight"]   = df.BMXWT
df["smoke_ever"] = yn(df.SMQ020)
# SMQ040 is only asked of ever-smokers; never-smokers are current non-smokers by definition
df["smoke_now"]  = df.SMQ040.map({1: 1, 2: 1, 3: 0})
df.loc[df.smoke_ever == 0, "smoke_now"] = 0
df["a1c"]      = df.LBXGH
df["dx_diab"]  = df.DIQ010.map({1: 1, 2: 0, 3: 0})   # 3 = borderline, treated as not diagnosed

for src, name in [("MCQ160A", "arthritis"), ("MCQ160B", "chf"), ("MCQ160C", "chd"),
                  ("MCQ160E", "mi"), ("MCQ160F", "stroke"), ("MCQ220", "cancer")]:
    df[name] = yn(df[src]) if src in df.columns else np.nan

# ---------- analytic sample ----------
df["adult"]    = df.age >= 20
df["examined"] = (df.RIDSTATR == 2) & df.WTMEC2YR.gt(0)
df["pregnant"] = (df.RIDEXPRG == 1)
ana = df[df.adult & df.examined & ~df.pregnant].copy()

# ---------- candidate outcomes ----------
ana["y_undiag_diabetes"] = np.where(ana.dx_diab == 0, (ana.a1c >= 6.5).astype(float), np.nan)
ana["y_obesity"]         = (ana.bmi >= 30).astype(float).where(ana.bmi.notna())
ana["y_central_obesity"] = np.where(ana.waist.notna() & ana.female.notna(),
                                    ((ana.female == 1) & (ana.waist > 88)) |
                                    ((ana.female == 0) & (ana.waist > 102)), np.nan)
ana["y_htn_measured"]    = np.where(ana.SBP.notna() & ana.DBP.notna(),
                                    ((ana.SBP >= 130) | (ana.DBP >= 80)).astype(float), np.nan)
for n in ["arthritis", "chf", "chd", "mi", "stroke", "cancer"]:
    ana[f"y_{n}"] = ana[n]

cols = (["SEQN", "cycle", "midyear", "bp_method", "SDMVPSU", "SDMVSTRA", "WTMEC2YR",
         "age", "female", "race", "educ_hs", "pir", "bmi", "waist", "height", "weight",
         "smoke_ever", "smoke_now", "a1c", "dx_diab", "SBP", "DBP"]
        + [c for c in ana.columns if c.startswith("y_")])
ana = ana[cols]
ana.to_csv(f"{OUT}/analytic.csv.gz", index=False, compression="gzip")

# ---------- structural report only (no predictor-outcome associations) ----------
print(f"rows: {len(ana):,}   cycles: {ana.cycle.nunique()}")
print()
print("per cycle: n, and non-missing count for each candidate outcome")
ycols = [c for c in ana.columns if c.startswith("y_")]
rep = ana.groupby("cycle").agg(n=("SEQN", "size"), **{c: (c, "count") for c in ycols})
print(rep.to_string())
print()
print("predictor missingness overall (%)")
pred = ["age", "female", "race", "educ_hs", "pir", "bmi", "waist",
        "smoke_ever", "smoke_now", "SBP"]
print((ana[pred].isna().mean() * 100).round(1).to_string())
