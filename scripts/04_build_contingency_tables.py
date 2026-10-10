"""
04_build_contingency_tables.py

Builds the 2x2 counts that script 05 turns into ROR and Information Component estimates.

Usage (from the repository root, after scripts 01, 02, 02b, 03 and 03b):
    python3 scripts/04_build_contingency_tables.py validation
    python3 scripts/04_build_contingency_tables.py extension

Drug groups, from the fluoroquinolone cases of script 02 (sole-suspect rule A or B, duplicate-removal
variant, and with or without topical-only cases): all fluoroquinolones, each of the six ingredients,
and combination (two or more fluoroquinolones).

Events: any PN24 term, any PN7 term, and each of the 24 terms individually (script 03).

Backgrounds (the comparison population; the paper only says "the full FAERS dataset"):
  all              every retained report in the pull (includes the fluoroquinolone reports)
  excl_fq_any      every retained report with no fluoroquinolone in any role
  no_concomitant   every retained report with no drug row in the concomitant role (C)
  paper_implied    all reports, but with the event fraction fixed at 1.826%, the value implied by
                   Table 1 of Shamim et al. 2026 (PN24 only)
The backgrounds 'all' and 'paper_implied' contain the drug group. 'no_concomitant' contains it under
rule B (the group has no other drug at all, apart from about 100 cases with a second fluoroquinolone
in a non-suspect role) but not under rule A. 'excl_fq_any' never contains it. Where the group is not
contained it is added back, so the comparison population is background + group.
Individual terms use the 'all' background only.

Columns of contingency.csv:
  O            reports with the drug group and the event
  n_drug       reports with the drug group
  n_event      reports with the event in the comparison population
  N            reports in the comparison population
  a, b, c, d   2x2 cells: a = O, b = n_drug - O, c = n_event - O, d = N - a - b - c

Output in data/interim/<pull>/:
  backgrounds.csv    reports and event reports for each background
  contingency.csv    one row per pull, rule, duplicate variant, topical scenario, group, event, background
"""
import importlib.util, sys
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PAPER_FRACTION = 0.01826


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s03 = load("s03", "03_filter_neuropathy_reactions.py")
ING, PN24, PN7 = s03.ING, s03.PN24, s03.PN7


def backgrounds(d, demo, fq):
    """Reports and PN events for each background, counted quarter by quarter."""
    fq_by_q = {q: set(g["rid"]) for q, g in fq.groupby("quarter")}
    acc = {k: [0, 0, 0] for k in ("all", "excl_fq_any", "no_concomitant")}
    s24, s7 = set(PN24), set(PN7)
    for q, g in demo.groupby("quarter"):
        rids = set(g["rid"])
        dr = pd.read_parquet(d / "drug" / f"{q}.parquet", columns=["rid", "role_cod"])
        conc = set(dr.loc[dr["role_cod"] == "C", "rid"])
        r = pd.read_parquet(d / "reac" / f"{q}.parquet", columns=["rid", "pt"]).drop_duplicates()
        r["pn"] = s03.normalise(r["pt"])
        h24 = set(r.loc[r["pn"].isin(s24), "rid"]) & rids
        h7 = set(r.loc[r["pn"].isin(s7), "rid"]) & rids
        fqs = fq_by_q.get(q, set())
        for name, pop in (("all", rids), ("excl_fq_any", rids - fqs), ("no_concomitant", rids - conc)):
            acc[name][0] += len(pop)
            acc[name][1] += len(pop & h24)
            acc[name][2] += len(pop & h7)
        print(f"{q}: {len(rids):,} reports", flush=True)
    bg = pd.DataFrame([(k, *v) for k, v in acc.items()], columns=["background", "N", "n_pn24", "n_pn7"])
    n = int(bg.loc[bg["background"] == "all", "N"].iloc[0])
    bg.loc[len(bg)] = ["paper_implied", n, round(PAPER_FRACTION * n), float("nan")]
    return bg


def row(base, O, n_drug, n_event, N, contains_group):
    if not contains_group:
        N, n_event = N + n_drug, n_event + O
    return {**base, "O": O, "n_drug": n_drug, "n_event": n_event, "N": N, "a": O, "b": n_drug - O,
            "c": n_event - O, "d": N - n_drug - (n_event - O)}


def main(pull):
    d = ROOT / "data" / "interim" / pull
    for f in ("demo.parquet", "fq_cases.parquet", "fq_cases_pn.parquet", "fq_topical.parquet",
              "fq_pn_terms.parquet", "background_counts.csv"):
        if not (d / f).exists():
            sys.exit(f"Missing {d / f}. Run scripts 01, 02, 02b, 03 and 03b first.")
    demo = pd.read_parquet(d / "demo.parquet", columns=["rid", "quarter"])
    fq = pd.read_parquet(d / "fq_cases.parquet", columns=["rid", "quarter"])
    bg = backgrounds(d, demo, fq)
    bg.to_csv(d / "backgrounds.csv", index=False)
    term_bg = pd.read_csv(d / "background_counts.csv").set_index("term")["reports"]
    N_all = int(bg.loc[bg["background"] == "all", "N"].iloc[0])

    t = pd.read_parquet(d / "fq_cases_pn.parquet").merge(
        pd.read_parquet(d / "fq_topical.parquet")[["rid", "quarter", "topical_only_any"]],
        on=["rid", "quarter"], how="left")
    t["topical_only_any"] = t["topical_only_any"].fillna(False).astype(bool)
    hits = pd.read_parquet(d / "fq_pn_terms.parquet").assign(v=True)
    wide = hits.set_index(["rid", "quarter", "pt"])["v"].unstack("pt", fill_value=False).reindex(
        columns=PN24, fill_value=False).astype(bool)
    t = t.merge(wide, left_on=["rid", "quarter"], right_index=True, how="left")
    t[PN24] = t[PN24].fillna(False).astype(bool)

    n_fq = sum(t[f"{i}_psss"].astype(int) for i in ING)
    out = []
    for rule, sole in (("A_no_other_suspect", ~t["has_nonfq_psss"]), ("B_no_other_drug_any_role", ~t["has_nonfq_any"])):
        for var, keep in (("missing_as_value", t["keep_missing_as_value"]), ("conservative", t["keep_conservative"])):
            for scen, sm in (("all_cases", pd.Series(True, index=t.index)), ("excl_topical_only", ~t["topical_only_any"])):
                m = (n_fq > 0) & sole & keep & sm
                groups = {"all": m, "combination": m & (n_fq >= 2)}
                groups.update({i: m & (n_fq == 1) & t[f"{i}_psss"] for i in ING})
                for g, k in groups.items():
                    base = {"pull": pull, "rule": rule, "dedup": var, "topical": scen, "group": g}
                    nd = int(k.sum())
                    for _, b in bg.iterrows():
                        for ev, col, col_bg in (("pn24", "pn24", "n_pn24"), ("pn7", "pn7", "n_pn7")):
                            if pd.isna(b[col_bg]):
                                continue
                            out.append(row({**base, "event": ev, "background": b["background"]},
                                           int((k & t[col]).sum()), nd, int(b[col_bg]), int(b["N"]),
                                           b["background"] in ("all", "paper_implied")
                                           or (b["background"] == "no_concomitant" and rule.startswith("B"))))
                    for term in PN24:
                        out.append(row({**base, "event": term, "background": "all"}, int((k & t[term]).sum()), nd,
                                       int(term_bg[term]), N_all, True))
    res = pd.DataFrame(out)
    res.to_csv(d / "contingency.csv", index=False)

    bad = int(((res[["a", "b", "c", "d"]] < 0).any(axis=1)).sum())
    if bad:
        print(f"WARNING: {bad:,} rows have a negative cell; check the background definitions.")
    print("\nBackgrounds:")
    print(bg.to_string(index=False))
    sel = res[(res["rule"] == "B_no_other_drug_any_role") & (res["dedup"] == "missing_as_value") &
              (res["topical"] == "all_cases") & (res["event"] == "pn24") & (res["background"] == "all")]
    print("\nRule B, duplicate variant 1, all cases, PN24, background 'all':")
    print(sel[["group", "O", "n_drug", "n_event", "N"]].to_string(index=False))
    print(f"\n{len(res):,} rows written to {d / 'contingency.csv'}")


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in s03.PULLS:
        sys.exit("usage: 04_build_contingency_tables.py [validation|extension]")
    main(sys.argv[1])
