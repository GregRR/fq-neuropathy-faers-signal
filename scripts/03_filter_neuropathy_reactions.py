"""
03_filter_neuropathy_reactions.py

Flags peripheral neuropathy (PN) related reactions on the fluoroquinolone cases and counts the same
terms across every report in the pull (the background used for the disproportionality step).

Usage (from the repository root, after scripts 01, 02 and 02b):
    python3 scripts/03_filter_neuropathy_reactions.py validation
    python3 scripts/03_filter_neuropathy_reactions.py extension

Two reaction sets, both MedDRA preferred terms as printed in Figure 2 of Shamim et al. 2026:
  PN24  the 24 PN-related terms shown in the paper's Figure 2 (terms with at least three reports).
        The paper selected 28 PN-related terms (its Table S9). Four terms had fewer than three
        reports and are not named in the main text, so they are missing here. Table 1 of the paper
        counts cases with any of the 28, so PN24 can only approximate that table.
  PN7   the seven terms the paper retained as FQ-associated signals (its Figure 2D).
A report with the same term spelled with US English (hypoesthesia, small fiber neuropathy) is
matched after a spelling normalisation; every raw spelling matched this way is listed in the audit.

Output in data/interim/<pull>/:
  fq_cases_pn.parquet     one row per fluoroquinolone case: the script 02 flags, the script 02b keep
                          flags, pn24, pn7
  fq_pn_terms.parquet     rid, quarter, pt (normalised) for each matched term on a fluoroquinolone case
  background_counts.csv   reports in the pull, reports with each term, reports with any PN24 / PN7 term
  pn_spelling_variants.csv  raw spellings matched only after normalisation
  pn_near_miss.csv        other terms that look related (not in either set), by report count, for review
  03_summary.json         case counts per ingredient with and without PN24, under both sole-suspect
                          readings and both duplicate-removal variants, next to Table 1 of the paper
"""
import json, re, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PULLS = ("validation", "extension")
ING = ["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin", "gemifloxacin", "delafloxacin"]

PN7 = ["paraesthesia", "neuropathy peripheral", "hypoaesthesia", "burning sensation",
       "sensory disturbance", "small fibre neuropathy", "hyperaesthesia"]
PN24 = PN7 + [
    "polyneuropathy", "sensory loss", "dysaesthesia", "unresponsive to stimuli",
    "electric shock sensation", "peripheral sensory neuropathy", "peripheral sensorimotor neuropathy",
    "allodynia", "sciatic nerve neuropathy", "mononeuropathy", "slow response to stimuli",
    "hyperpathia", "acute polyneuropathy", "toxic neuropathy", "demyelinating polyneuropathy",
    "peripheral motor neuropathy", "axonal neuropathy"]
assert len(PN7) == 7 and len(PN24) == 24 and len(set(PN24)) == 24
NEAR = re.compile(r"aesthesia|neuropath|burning|sensory|sensation|numb|tingl|pins and needles|nerve|"
                  r"paralysis|paresis|hyperpathia|hypoalgesia|neuralgia|neuritis|dysesthesia|neurotox")

# Table 1 of Shamim et al. 2026: cases with PN-related AEs / total cases
PAPER_PN = {"all": (5588, 37778), "ciprofloxacin": (2291, 14434), "levofloxacin": (1942, 13598),
            "moxifloxacin": (642, 7439), "ofloxacin": (46, 588), "gemifloxacin": (11, 142),
            "delafloxacin": (4, 72), "combination": (652, 1505)}
# Figure 2A / Section 3.2 of the paper: cases with the term, all fluoroquinolones
PAPER_TERMS = {"paraesthesia": 2095, "neuropathy peripheral": 1986, "hypoaesthesia": 1637,
               "burning sensation": 1066}


def normalise(s):
    """Lower case, trim, UK spelling. Done on the unique values for speed."""
    u = pd.Series(s.unique())
    n = (u.str.strip().str.lower().str.replace("fiber", "fibre", regex=False)
          .str.replace(r"(?<!a)esthesia", "aesthesia", regex=True))
    return s.map(dict(zip(u, n)))


def summarize(t):
    n_fq = sum(t[f"{i}_psss"].astype(int) for i in ING)
    out = {}
    for rule, sole in (("A_no_other_suspect", ~t["has_nonfq_psss"]), ("B_no_other_drug_any_role", ~t["has_nonfq_any"])):
        for var, keep in (("missing_as_value", t["keep_missing_as_value"]), ("conservative", t["keep_conservative"])):
            m = (n_fq > 0) & sole & keep
            groups = {"all": m, "combination": m & (n_fq >= 2)}
            groups.update({i: m & (n_fq == 1) & t[f"{i}_psss"] for i in ING})
            res = {g: {"cases": int(k.sum()), "pn24": int((k & t["pn24"]).sum()), "pn7": int((k & t["pn7"]).sum())}
                   for g, k in groups.items()}
            res["terms_all_fq"] = {p: int((m & t[f"term_{j}"]).sum()) for j, p in enumerate(PN24)}
            out[f"{rule}__{var}"] = res
    return out


def main(pull):
    d = ROOT / "data" / "interim" / pull
    need = [d / "fq_cases.parquet", d / "fq_cases_dedup.parquet", d / "demo.parquet"]
    for f in need:
        if not f.exists():
            sys.exit(f"Missing {f}. Run scripts 01, 02 and 02b first.")
    fq = pd.read_parquet(d / "fq_cases.parquet")
    dd = pd.read_parquet(d / "fq_cases_dedup.parquet", columns=["rid", "quarter", "caseid", "keep_missing_as_value",
                                                                "keep_conservative", "group_size"])
    t = fq.merge(dd, on=["rid", "quarter"], how="left", validate="one_to_one")
    if t["keep_missing_as_value"].isna().any():
        sys.exit("Some fluoroquinolone cases have no duplicate-removal flag. Re-run scripts 02 and 02b.")
    n_reports = len(pd.read_parquet(d / "demo.parquet", columns=["rid"]))

    s24, s7 = set(PN24), set(PN7)
    fq_by_q = {q: set(g["rid"]) for q, g in t.groupby("quarter")}
    term_counts = pd.Series(0, index=PN24, dtype="int64")
    any24 = any7 = 0
    variants, near, fq_hits = [], [], []
    for p in sorted((d / "reac").glob("*.parquet")):
        q = p.stem
        r = pd.read_parquet(p, columns=["rid", "pt"]).drop_duplicates()
        r["pn"] = normalise(r["pt"])
        hit = r[r["pn"].isin(s24)].drop_duplicates(["rid", "pn"])
        term_counts = term_counts.add(hit.groupby("pn").size(), fill_value=0)
        any24 += hit["rid"].nunique()
        any7 += hit.loc[hit["pn"].isin(s7), "rid"].nunique()
        v = hit[hit["pt"].str.strip().str.lower() != hit["pn"]]
        if len(v):
            variants.append(v.groupby(["pt", "pn"]).size().rename("reports").reset_index())
        miss = r[~r["pn"].isin(s24) & r["pn"].str.contains(NEAR)]
        near.append(miss.groupby("pn").size().rename("reports").reset_index())
        keep = fq_by_q.get(q)
        if keep:
            fh = hit[hit["rid"].isin(keep)][["rid", "pn"]].rename(columns={"pn": "pt"})
            fq_hits.append(fh.assign(quarter=q))
        print(f"{q}: {r['rid'].nunique():,} reports with reactions", flush=True)

    hits = pd.concat(fq_hits, ignore_index=True)[["rid", "quarter", "pt"]]
    hits.to_parquet(d / "fq_pn_terms.parquet", index=False)
    key = hits.assign(v=True).set_index(["rid", "quarter", "pt"])["v"]
    wide = key.unstack("pt", fill_value=False).reindex(columns=PN24, fill_value=False).astype(bool)
    wide.columns = [f"term_{j}" for j in range(len(PN24))]
    t = t.merge(wide, left_on=["rid", "quarter"], right_index=True, how="left")
    tc = [f"term_{j}" for j in range(len(PN24))]
    t[tc] = t[tc].fillna(False).astype(bool)
    t["pn24"] = t[tc].any(axis=1)
    t["pn7"] = t[[f"term_{j}" for j in range(7)]].any(axis=1)
    t[[c for c in t.columns if c not in tc]].to_parquet(d / "fq_cases_pn.parquet", index=False)

    bg = pd.DataFrame({"term": ["__reports_in_pull__", "__any_PN24__", "__any_PN7__"] + PN24,
                       "reports": [n_reports, any24, any7] + [int(term_counts[x]) for x in PN24]})
    bg.to_csv(d / "background_counts.csv", index=False)
    pd.concat(variants).groupby(["pt", "pn"], as_index=False)["reports"].sum().to_csv(
        d / "pn_spelling_variants.csv", index=False) if variants else pd.DataFrame(columns=["pt", "pn", "reports"]).to_csv(
        d / "pn_spelling_variants.csv", index=False)
    pd.concat(near).groupby("pn", as_index=False)["reports"].sum().sort_values("reports", ascending=False).to_csv(
        d / "pn_near_miss.csv", index=False)

    res = summarize(t)
    (d / "03_summary.json").write_text(json.dumps({"reports_in_pull": n_reports, "pn24_reports": any24,
                                                  "pn7_reports": any7, "paper_table1": PAPER_PN, **res}, indent=2))

    cols = ["all"] + ING + ["combination"]
    print(f"\nReports in pull: {n_reports:,}; with any PN24 term: {any24:,}; any PN7 term: {any7:,}")
    print("\nCases with a PN24 term / cases (compare Table 1):")
    for k, v in res.items():
        print(f"  {k:46s}: " + ", ".join(f"{c}={v[c]['pn24']:,}/{v[c]['cases']:,}" for c in cols))
    print("  " + "paper".ljust(46) + ": " + ", ".join(f"{c}={PAPER_PN[c][0]:,}/{PAPER_PN[c][1]:,}" for c in cols))
    k = "B_no_other_drug_any_role__missing_as_value"
    print(f"\nCases with each of the four most common terms ({k}, all fluoroquinolones):")
    for term, n in PAPER_TERMS.items():
        print(f"  {term:24s} {res[k]['terms_all_fq'][term]:>6,}  (paper {n:,})")
    print(f"\nWritten to {d}")


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in PULLS:
        sys.exit("usage: 03_filter_neuropathy_reactions.py [validation|extension]")
    main(sys.argv[1])
