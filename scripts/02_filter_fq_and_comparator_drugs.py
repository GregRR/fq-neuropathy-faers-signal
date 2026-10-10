"""
02_filter_fq_and_comparator_drugs.py

Flags the six FDA-approved fluoroquinolones (ciprofloxacin, levofloxacin, moxifloxacin,
ofloxacin, gemifloxacin, delafloxacin) in the DRUG rows written by script 01 and builds
one row per case that reports at least one of them, in any role.

Usage (from the repository root):
    python3 scripts/02_filter_fq_and_comparator_drugs.py validation
    python3 scripts/02_filter_fq_and_comparator_drugs.py extension

Drug identification: DRUGNAME and PROD_AI (active ingredient, available from 2014q3) are
searched for generic-name stems and for the U.S. brand names Cipro, Proquin, Levaquin,
Avelox, Floxin, Factive and Baxdela (brand list: drugs.com quinolone class page).
Legacy files have no PROD_AI, so those rows rely on DRUGNAME alone. To close that gap, a
lookup is built from rows where the FDA's own PROD_AI identifies an ingredient but DRUGNAME
matched none of the patterns (for example CRAVIT, TAVANIC, CIPROXIN). A name enters the
lookup only if it appears at least MIN_ROWS times and PROD_AI maps it to one ingredient.
The lookup is applied to the DRUGNAME of every row, and is saved for review as
fq_name_map.csv.
Audit tables are written next to the output so name coverage can be reviewed:
  fq_name_map.csv             the PROD_AI-derived DRUGNAME lookup described above
  fq_name_audit.csv           every DRUGNAME matched to each ingredient, with row counts
  fq_via_prod_ai_only.csv     names matched only through PROD_AI (candidate extra patterns)
  quinolone_stem_unmatched.csv names containing FLOX/QUINOLON not matched to the six
  fq_route_audit.csv          ROUTE values per ingredient (topical vs systemic review)

Output in data/interim/<pull>/:
  fq_cases.parquet   rid, <ingredient>_psss and <ingredient>_any (bool, role PS/SS or any
                     role), has_nonfq_psss, has_nonfq_any, has_other_quinolone
  02_summary.json    case counts under two readings of "sole suspected drug":
      A: at least one FQ as PS/SS and no non-FQ drug as PS/SS (concomitant drugs allowed)
      B: at least one FQ as PS/SS and no non-FQ drug of any role
    each split into single FQ (exactly one FQ ingredient as PS/SS) and combination (two or
    more), for comparison with Table 1 of Shamim et al. 2026.
"""
import json, re, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
FQ = {
    "ciprofloxacin": r"\bCIPROFLOX|\bCIPRO\b|\bPROQUIN\b",
    "levofloxacin": r"\bLEVOFLOX|\bLEVAQUIN\b",
    "moxifloxacin": r"\bMOXIFLOX|\bAVELOX\b",
    # \b before OFLOXACIN so that LEVOFLOXACIN and CIPROFLOXACIN do not match it
    "ofloxacin": r"\bOFLOXACIN|\bFLOXIN\b",
    "gemifloxacin": r"\bGEMIFLOX|\bFACTIVE\b",
    "delafloxacin": r"\bDELAFLOX|\bBAXDELA\b",
}
MIN_ROWS = 3
STEM = re.compile(r"FLOXACIN|FLOXAC|QUINOLON")
RX = {k: re.compile(v) for k, v in FQ.items()}
PAPER = {"all": 37778, "ciprofloxacin": 14434, "levofloxacin": 13598, "moxifloxacin": 7439,
         "ofloxacin": 588, "gemifloxacin": 142, "delafloxacin": 72, "combination": 1505}


def classify_names(uniq, name_map=None):
    """uniq: DataFrame of distinct (drugname, prod_ai). Returns bool columns per ingredient."""
    out = pd.DataFrame(index=uniq.index)
    for ing, rx in RX.items():
        out[ing] = uniq["drugname"].str.contains(rx) | uniq["prod_ai"].str.contains(rx)
        if name_map:
            out[ing] |= uniq["drugname"].isin(name_map.get(ing, ()))
    return out


def build_name_map(parts):
    """Names the FDA's PROD_AI links to exactly one ingredient but the patterns do not match."""
    frames = []
    for p in parts:
        df = pd.read_parquet(p, columns=["drugname", "prod_ai"]).fillna("")
        df = df[(df["prod_ai"] != "") & (df["drugname"] != "")]
        g = df.groupby(["drugname", "prod_ai"]).size().rename("rows").reset_index()
        flags = classify_names(g[["drugname", "prod_ai"]])
        for ing, rx in RX.items():
            by_ai = g["prod_ai"].str.contains(rx)
            by_name = g["drugname"].str.contains(rx)
            sel = g[by_ai & ~by_name]
            if len(sel):
                frames.append(sel.assign(ingredient=ing)[["drugname", "ingredient", "rows"]])
    if not frames:
        return {}, pd.DataFrame(columns=["drugname", "ingredient", "rows"])
    tab = pd.concat(frames).groupby(["drugname", "ingredient"], as_index=False)["rows"].sum()
    n_ing = tab.groupby("drugname")["ingredient"].transform("nunique")
    keep = tab[(tab["rows"] >= MIN_ROWS) & (n_ing == 1)].sort_values("rows", ascending=False)
    name_map = {ing: set(g["drugname"]) for ing, g in keep.groupby("ingredient")}
    return name_map, keep


def process_partition(path, name_map, name_counts, via_ai, stem_counts, route_counts):
    df = pd.read_parquet(path)
    df["drugname"] = df["drugname"].fillna("")
    df["prod_ai"] = df["prod_ai"].fillna("")
    uniq = df[["drugname", "prod_ai"]].drop_duplicates().reset_index(drop=True)
    flags = classify_names(uniq, name_map)
    any_fq = flags.any(axis=1)
    uniq = pd.concat([uniq, flags], axis=1)
    uniq["any_fq"] = any_fq
    uniq["stem_unmatched"] = (uniq["drugname"].str.contains(STEM) | uniq["prod_ai"].str.contains(STEM)) & ~any_fq
    df = df.merge(uniq, on=["drugname", "prod_ai"], how="left")

    ps = df["role_cod"].isin(["PS", "SS"])
    ings = list(FQ)
    for ing in ings:
        df[f"{ing}_any"] = df[ing]
        df[f"{ing}_psss"] = df[ing] & ps
    df["has_nonfq_any"] = ~df["any_fq"]
    df["has_nonfq_psss"] = ~df["any_fq"] & ps
    df["has_other_quinolone"] = df["stem_unmatched"]
    cols = [f"{i}_{s}" for i in ings for s in ("any", "psss")] + ["has_nonfq_any", "has_nonfq_psss", "has_other_quinolone"]
    cases = df.groupby("rid", sort=False)[cols].max()
    cases = cases[cases[[f"{i}_any" for i in ings]].any(axis=1)].reset_index()

    # audit tables
    m = df[df["any_fq"]]
    for ing in ings:
        sub = m[m[ing]]
        name_counts.append(sub.groupby("drugname").size().rename("rows").reset_index().assign(ingredient=ing))
        route_counts.append(sub.groupby("route").size().rename("rows").reset_index().assign(ingredient=ing))
    ai_only = m[~m["drugname"].str.contains("|".join(FQ.values())) & (m["prod_ai"] != "")]
    via_ai.append(ai_only.groupby(["drugname", "prod_ai"]).size().rename("rows").reset_index())
    stem_counts.append(df[df["stem_unmatched"]].groupby("drugname").size().rename("rows").reset_index())
    return cases


def summarize(fq):
    ings = list(FQ)
    n_fq = sum(fq[f"{i}_psss"].astype(int) for i in ings)
    has = n_fq > 0
    res = {}
    for label, sole in (("A_no_other_suspect", ~fq["has_nonfq_psss"]), ("B_no_other_drug_any_role", ~fq["has_nonfq_any"])):
        keep = has & sole
        r = {"all": int(keep.sum()), "combination": int((keep & (n_fq >= 2)).sum())}
        for i in ings:
            r[i] = int((keep & (n_fq == 1) & fq[f"{i}_psss"]).sum())
        r["paper_table1"] = PAPER
        res[label] = r
    res["fq_suspect_any_other_drugs"] = {"all": int(has.sum()),
                                           **{i: int((has & (n_fq == 1) & fq[f"{i}_psss"]).sum()) for i in ings},
                                           "combination": int((has & (n_fq >= 2)).sum())}
    res["fq_cases_any_role"] = int(len(fq))
    res["cases_with_other_quinolone_name"] = int(fq["has_other_quinolone"].sum())
    return res


def main(pull):
    d = ROOT / "data" / "interim" / pull
    parts = sorted((d / "drug").glob("*.parquet"))
    if not parts:
        sys.exit(f"No drug partitions in {d / 'drug'}. Run script 01 first.")
    name_map, name_map_tab = build_name_map(parts)
    name_map_tab.to_csv(d / "fq_name_map.csv", index=False)
    print(f"Name lookup from PROD_AI: {len(name_map_tab):,} names", flush=True)
    name_counts, via_ai, stem_counts, route_counts, cases = [], [], [], [], []
    for p in parts:
        cases.append(process_partition(p, name_map, name_counts, via_ai, stem_counts, route_counts))
        print(f"{p.stem}: {len(cases[-1]):,} cases with a fluoroquinolone", flush=True)
    fq = pd.concat(cases, ignore_index=True)
    fq.to_parquet(d / "fq_cases.parquet", index=False)

    def dump(frames, by, name):
        out = pd.concat(frames).groupby(by, as_index=False)["rows"].sum().sort_values("rows", ascending=False)
        out.to_csv(d / name, index=False)
    dump(name_counts, ["ingredient", "drugname"], "fq_name_audit.csv")
    dump(via_ai, ["drugname", "prod_ai"], "fq_via_prod_ai_only.csv")
    dump(stem_counts, ["drugname"], "quinolone_stem_unmatched.csv")
    dump(route_counts, ["ingredient", "route"], "fq_route_audit.csv")

    res = summarize(fq)
    res["name_lookup_names"] = int(len(name_map_tab))
    (d / "02_summary.json").write_text(json.dumps(res, indent=2))
    print("\nCase counts (single FQ by ingredient, then combination):")
    for k in ("A_no_other_suspect", "B_no_other_drug_any_role"):
        print(f"  {k}: " + ", ".join(f"{a}={res[k][a]:,}" for a in ("all", *FQ, "combination")))
    print("  paper : " + ", ".join(f"{a}={PAPER[a]:,}" for a in ("all", *FQ, "combination")))
    print(f"\nSummary and audit tables written to {d}")


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("validation", "extension"):
        sys.exit("usage: 02_filter_fq_and_comparator_drugs.py [validation|extension]")
    main(sys.argv[1])
