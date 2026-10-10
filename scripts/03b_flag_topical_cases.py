"""
03b_flag_topical_cases.py

Sensitivity step: flags fluoroquinolone cases whose suspect fluoroquinolone rows are all topical
(eye, ear or skin products), and shows how excluding them changes the case counts.

The paper does not say whether topical products were included. Topical products rarely cause
peripheral neuropathy, so including them adds cases to the denominators without adding neuropathy cases.

Usage (from the repository root, after scripts 01, 02, 02b and 03):
    python3 scripts/03b_flag_topical_cases.py validation
    python3 scripts/03b_flag_topical_cases.py extension

A fluoroquinolone drug row (role PS or SS) is topical when
  by route:  the route is ophthalmic, intraocular, intracorneal, (sub)conjunctival, intra-aural,
             auricular (otic), topical, cutaneous or occlusive dressing technique
  by name:   the drug name contains OPHTH, OTIC, EYE, EAR, DROP, AUGEN, COLLYR, AURICUL, OINTMENT,
             UNGUENTO or CREAM, or is one of the brands listed as ophthalmic or otic fluoroquinolones in
             the Missouri HealthNet ophthalmic fluoroquinolone review and the Nevada Medicaid otic
             fluoroquinolone review (Vigamox, Moxeza, Ciloxan, Ocuflox, Cetraxal, Otiprio, Ciprodex,
             Otovel, Cipro HC)
A case is topical-only when every one of its suspect fluoroquinolone rows is topical. Rows with a blank
or unknown route and no topical name count as not topical, so the flag is conservative.

Output in data/interim/<pull>/:
  fq_topical.parquet       rid, quarter, suspect fluoroquinolone rows, topical rows (route, route or name),
                           topical_only_route, topical_only_any
  topical_name_audit.csv   drug names behind the flags, with route and row counts, for review
  03b_summary.json         case counts with and without the topical-only cases, next to Table 1
"""
import importlib.util, json, re, sys
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s02 = load("s02", "02_filter_fq_and_comparator_drugs.py")
s03 = load("s03", "03_filter_neuropathy_reactions.py")
ING = s03.ING
TOPICAL_ROUTE = {"OPHTHALMIC", "INTRAOCULAR", "INTRACORNEAL", "SUNCONJUNCTIVAL", "SUBCONJUNCTIVAL",
                 "INTRA-AURAL", "AURICULAR (OTIC)", "TOPICAL", "CUTANEOUS", "OCCLUSIVE DRESSING TECHNIQUE"}
TOPICAL_NAME = re.compile(r"OPHTH|\bOTIC\b|\bEYES?\b|\bEARS?\b|DROP|AUGEN|COLLYR|AURICUL|OINTMENT|UNGUENTO|CREAM|"
                          r"\bVIGAMOX|\bMOXEZA|\bCILOXAN|\bOCUFLOX|\bCETRAXAL|\bOTIPRIO|\bCIPRODEX|\bOTOVEL|\bCIPRO HC\b")


def summarize(t, mask):
    n_fq = sum(t[f"{i}_psss"].astype(int) for i in ING)
    out = {}
    for rule, sole in (("A_no_other_suspect", ~t["has_nonfq_psss"]), ("B_no_other_drug_any_role", ~t["has_nonfq_any"])):
        for var, keep in (("missing_as_value", t["keep_missing_as_value"]), ("conservative", t["keep_conservative"])):
            m = (n_fq > 0) & sole & keep & mask
            groups = {"all": m, "combination": m & (n_fq >= 2)}
            groups.update({i: m & (n_fq == 1) & t[f"{i}_psss"] for i in ING})
            out[f"{rule}__{var}"] = {g: {"cases": int(k.sum()), "pn24": int((k & t["pn24"]).sum())}
                                     for g, k in groups.items()}
    return out


def main(pull):
    d = ROOT / "data" / "interim" / pull
    if not (d / "fq_cases_pn.parquet").exists():
        sys.exit(f"Missing {d / 'fq_cases_pn.parquet'}. Run scripts 01, 02, 02b and 03 first.")
    lookup_parts = sorted(q for pl in ("validation", "extension")
                          for q in (ROOT / "data" / "interim" / pl / "drug").glob("*.parquet"))
    name_map, _ = s02.build_name_map(lookup_parts)
    rows, audit = [], []
    for p in sorted((d / "drug").glob("*.parquet")):
        df = pd.read_parquet(p, columns=["rid", "role_cod", "drugname", "prod_ai", "route"]).fillna("")
        df = df[df["role_cod"].isin(["PS", "SS"])]
        uniq = df[["drugname", "prod_ai"]].drop_duplicates().reset_index(drop=True)
        uniq["fq"] = s02.classify_names(uniq, name_map).any(axis=1)
        df = df.merge(uniq, on=["drugname", "prod_ai"], how="left")
        df = df[df["fq"]].copy()
        df["by_route"] = df["route"].isin(TOPICAL_ROUTE)
        df["by_name"] = df["drugname"].str.contains(TOPICAL_NAME)
        df["topical"] = df["by_route"] | df["by_name"]
        audit.append(df[df["topical"]].groupby(["drugname", "route", "by_route", "by_name"]).size()
                     .rename("rows").reset_index())
        g = df.groupby("rid").agg(fq_rows=("fq", "size"), topical_route_rows=("by_route", "sum"),
                                  topical_any_rows=("topical", "sum")).reset_index()
        rows.append(g.assign(quarter=p.stem))
        print(f"{p.stem}: {len(g):,} cases with a suspect fluoroquinolone row", flush=True)
    f = pd.concat(rows, ignore_index=True)
    f["topical_only_route"] = f["topical_route_rows"] == f["fq_rows"]
    f["topical_only_any"] = f["topical_any_rows"] == f["fq_rows"]
    f.to_parquet(d / "fq_topical.parquet", index=False)
    pd.concat(audit).groupby(["drugname", "route", "by_route", "by_name"], as_index=False)["rows"].sum() \
        .sort_values("rows", ascending=False).to_csv(d / "topical_name_audit.csv", index=False)

    t = pd.read_parquet(d / "fq_cases_pn.parquet").merge(
        f[["rid", "quarter", "topical_only_route", "topical_only_any"]], on=["rid", "quarter"], how="left")
    t[["topical_only_route", "topical_only_any"]] = t[["topical_only_route", "topical_only_any"]].fillna(False).astype(bool)
    scen = {"all_cases": pd.Series(True, index=t.index), "excl_topical_only_route": ~t["topical_only_route"],
            "excl_topical_only_route_or_name": ~t["topical_only_any"]}
    res = {k: summarize(t, m) for k, m in scen.items()}
    (d / "03b_summary.json").write_text(json.dumps({"paper_table1": s03.PAPER_PN, **res}, indent=2))

    cols = ["all"] + ING + ["combination"]
    print(f"\nTopical-only cases among all fluoroquinolone cases: by route {int(t['topical_only_route'].sum()):,}, "
          f"by route or name {int(t['topical_only_any'].sum()):,} of {len(t):,}")
    print("\nCases with a PN24 term / cases:")
    for var in ("missing_as_value", "conservative"):
        k = f"B_no_other_drug_any_role__{var}"
        for sc in scen:
            print(f"  B {var[:12]:12s} {sc:34s}: " + ", ".join(
                f"{c}={res[sc][k][c]['pn24']:,}/{res[sc][k][c]['cases']:,}" for c in cols))
    print("  " + "paper".ljust(48) + ": " + ", ".join(
        f"{c}={s03.PAPER_PN[c][0]:,}/{s03.PAPER_PN[c][1]:,}" for c in cols))
    print(f"\nWritten to {d}")


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in s03.PULLS:
        sys.exit("usage: 03b_flag_topical_cases.py [validation|extension]")
    main(sys.argv[1])
