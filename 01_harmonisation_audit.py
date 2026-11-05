"""
NHANES harmonisation audit, step 1.
Reads every staged .xpt, records row counts and variable inventories per cycle,
and reports which variables are present in which cycles.
No modelling, no outcome inspection.
"""
import glob, os, json, hashlib
import pandas as pd, pyreadstat

SRC = "/mnt/user-data/uploads/Downloads/NHANES_JAMIA/_frozen_v1"
OUT = "/home/claude/jamia2026/03_analysis/out"
os.makedirs(OUT, exist_ok=True)

CYCLE = {"": "1999-2000", "_B": "2001-2002", "_C": "2003-2004", "_D": "2005-2006",
         "_E": "2007-2008", "_F": "2009-2010", "_G": "2011-2012", "_H": "2013-2014",
         "_I": "2015-2016", "_J": "2017-2018", "_L": "2021-2023"}

def split_name(stem):
    # LAB10 / L10_B / L10_C are the 1999-2004 names for the lipid+glucose file
    for suf in ["_B","_C","_D","_E","_F","_G","_H","_I","_J","_L"]:
        if stem.endswith(suf):
            return stem[: -len(suf)], suf
    return stem, ""

rows, varmap = [], {}
for p in sorted(glob.glob(os.path.join(SRC, "*.xpt"))):
    stem = os.path.basename(p)[:-4]
    comp, suf = split_name(stem)
    cyc = CYCLE.get(suf)
    try:
        df, meta = pyreadstat.read_xport(p, encoding="latin1")
        cols = list(meta.column_names)
    except Exception:
        df = pd.read_sas(p, format="xport", encoding="latin1")
        cols = list(df.columns)
    sha = hashlib.sha256(open(p,"rb").read()).hexdigest()
    rows.append(dict(file=stem, component=comp, cycle=cyc, n_rows=len(df),
                     n_vars=len(cols), sha256=sha))
    for v in cols:
        varmap.setdefault((comp, v), set()).add(cyc)

inv = pd.DataFrame(rows).sort_values(["component","cycle"])
inv.to_csv(f"{OUT}/file_inventory.csv", index=False)

cycles_present = sorted({r["cycle"] for r in rows if r["cycle"]})
vm = []
for (comp, var), cycs in varmap.items():
    vm.append(dict(component=comp, variable=var, n_cycles=len(cycs),
                   cycles=";".join(sorted(cycs)),
                   in_all=len(cycs) == len(cycles_present)))
vmdf = pd.DataFrame(vm).sort_values(["component","n_cycles","variable"],
                                    ascending=[True,False,True])
vmdf.to_csv(f"{OUT}/variable_coverage.csv", index=False)

print("cycles present:", len(cycles_present))
for c in cycles_present: print("  ", c)
print()
print(inv.groupby("component")[["n_rows"]].agg(["count","min","max"]).to_string())
print()
print("variables present in ALL cycles, by component:")
for comp, g in vmdf[vmdf.in_all].groupby("component"):
    print(f"  {comp:6s} ({len(g):3d}): {', '.join(sorted(g.variable)[:28])}")
