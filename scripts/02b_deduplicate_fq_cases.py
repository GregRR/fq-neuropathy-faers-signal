"""
02b_deduplicate_fq_cases.py

Second-stage duplicate removal for fluoroquinolone cases, following the rule described by
Shamim et al. 2026 (Clin Transl Sci 19:e70541): reports with identical age, sex, country,
initial FDA receipt date, drug, indication and reaction are treated as duplicates and one is
kept. Script 01 already removed repeated versions of the same case (same CASEID); this step
catches the same event reported under different CASEIDs (for example by a manufacturer and
by a clinician).

Usage (from the repository root, after scripts 01 and 02):
    python3 scripts/02b_deduplicate_fq_cases.py validation
    python3 scripts/02b_deduplicate_fq_cases.py extension

Match fields (documented choices where the paper is not specific):
  age        AGE converted to years, rounded to 0.1
  sex        SEX / GNDR_COD
  country    occurrence country (OCCR_COUNTRY) when present, else reporter country
             (legacy files have only the reporter country)
  date       INIT_FDA_DT. Legacy reports have no INIT_FDA_DT; their FDA_DT is used instead.
  drug       the set of fluoroquinolone ingredients reported as primary or secondary suspect, plus the
             set of every other drug on the report in any role (PROD_AI when present, else the drug name)
  indication the set of INDI_PT terms across all drugs on the report (from the INDI files)
  reaction   the set of REAC PT terms
Two variants are written:
  keep_missing_as_value  a missing age, sex or country counts as a value, so two reports
                         that are both missing the same field can match
  keep_conservative      reports missing age, sex or country are never merged
Within a duplicate group the report with the lowest CASEID is kept.

Output in data/interim/<pull>/:
  fq_cases_dedup.parquet  rid, quarter, caseid, the two keep flags and the group size
  02b_summary.json        case counts by ingredient before and after, under both sole-suspect
                          readings from script 02, next to Table 1 of the paper
Duplicate removal is applied to fluoroquinolone cases only; the full-FAERS background used
later for the BCPNN is not deduplicated this way.
"""
import importlib.util, json, re, sys, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("s01", HERE / "01_parse_and_link_core_files.py")
s01 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s01)
spec2 = importlib.util.spec_from_file_location("s02", HERE / "02_filter_fq_and_comparator_drugs.py")
s02 = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(s02)

INDI_MEMBER = re.compile(r"^indi\d{2}q[1-4](?:_new)?\.txt$", re.I)
ING = list(s02.FQ)
PAPER = s02.PAPER


def read_indications(pull, rids_by_quarter):
    """INDI_PT terms per report for the given reports, read straight from the ZIPs."""
    start, end = s01.PULLS[pull]
    frames = []
    for y, q in s01.quarters(start, end):
        label = f"{y}q{q}"
        keep = rids_by_quarter.get(label)
        if not keep:
            continue
        schema = "L" if (y, q) <= (2012, 3) else "F"
        with zipfile.ZipFile(s01.zip_path(y, q)) as zf:
            names = [n for n in zf.namelist() if INDI_MEMBER.match(Path(n).name)]
            if len(names) != 1:
                sys.exit(f"{zf.filename}: expected one INDI file, found {names}")
            stats = {"lines_skipped": {}}
            df = s01.read_member(zf, names[0], stats, f"{label}:indi")
        idcol = "isr" if schema == "L" else "primaryid"
        if idcol not in df or "indi_pt" not in df:
            sys.exit(f"{label}: INDI columns {list(df.columns)} lack {idcol} / indi_pt")
        df["rid"] = schema + df[idcol]
        df = df[df["rid"].isin(keep)][["rid", "indi_pt"]].assign(quarter=label)
        df["indi_pt"] = df["indi_pt"].str.strip().str.lower()
        frames.append(df)
        print(f"indications {label}: {len(df):,} rows for {len(keep):,} cases", flush=True)
    return pd.concat(frames, ignore_index=True)


def other_drug_keys(pull, rids_by_quarter):
    """Set of non-fluoroquinolone drugs (any role) per report: PROD_AI when present, else the drug name."""
    lookup_parts = sorted(q for pl in ("validation", "extension")
                          for q in (ROOT_DIR / "data" / "interim" / pl / "drug").glob("*.parquet"))
    name_map, _ = s02.build_name_map(lookup_parts)
    frames = []
    for q, keep in rids_by_quarter.items():
        df = pd.read_parquet(ROOT_DIR / "data" / "interim" / pull / "drug" / f"{q}.parquet",
                             columns=["rid", "drugname", "prod_ai"]).fillna("")
        df = df[df["rid"].isin(keep)]
        uniq = df[["drugname", "prod_ai"]].drop_duplicates().reset_index(drop=True)
        uniq["fq"] = s02.classify_names(uniq, name_map).any(axis=1)
        uniq["ident"] = uniq["prod_ai"].where(uniq["prod_ai"] != "", uniq["drugname"])
        df = df.merge(uniq, on=["drugname", "prod_ai"], how="left")
        df = df[~df["fq"] & (df["ident"] != "")]
        frames.append(df[["rid", "ident"]].assign(quarter=q))
    d = pd.concat(frames, ignore_index=True)
    return d.groupby(["rid", "quarter"])["ident"].agg(lambda s: "|".join(sorted(set(s)))).rename("other_key").reset_index()


def set_key(df, col):
    df = df[df[col] != ""]
    return df.groupby(["rid", "quarter"])[col].agg(lambda s: "|".join(sorted(set(s)))).rename(col + "_key").reset_index()


def summarize(fq, keep):
    n_fq = sum(fq[f"{i}_psss"].astype(int) for i in ING)
    out = {}
    for label, sole in (("A_no_other_suspect", ~fq["has_nonfq_psss"]),
                        ("B_no_other_drug_any_role", ~fq["has_nonfq_any"])):
        for tag, k in (("before", pd.Series(True, index=fq.index)), ("after", keep)):
            m = (n_fq > 0) & sole & k
            r = {"all": int(m.sum()), "combination": int((m & (n_fq >= 2)).sum())}
            for i in ING:
                r[i] = int((m & (n_fq == 1) & fq[f"{i}_psss"]).sum())
            out[f"{label}_{tag}"] = r
    return out


def main(pull):
    d = ROOT_DIR / "data" / "interim" / pull
    for f in ("fq_cases.parquet", "demo.parquet"):
        if not (d / f).exists():
            sys.exit(f"Missing {d / f}. Run scripts 01 and 02 first.")
    fq = pd.read_parquet(d / "fq_cases.parquet")
    demo = pd.read_parquet(d / "demo.parquet")
    cols = ["rid", "caseid", "quarter", "schema", "age_years", "sex", "reporter_country",
            "occr_country", "init_fda_dt", "fda_dt"]
    if "quarter" not in fq:
        sys.exit("fq_cases.parquet has no quarter column. Re-run script 02.")
    c = fq.merge(demo[cols], on=["rid", "quarter"], how="left", validate="one_to_one")

    rids_by_q = {q: set(g) for q, g in c.groupby("quarter")["rid"]}
    reac_frames = []
    for q, keep in rids_by_q.items():
        r = pd.read_parquet(d / "reac" / f"{q}.parquet")
        reac_frames.append(r[r["rid"].isin(keep)].assign(quarter=q))
    reac = pd.concat(reac_frames, ignore_index=True)
    indi = read_indications(pull, rids_by_q)

    c = c.merge(set_key(reac, "pt"), on=["rid", "quarter"], how="left")
    c = c.merge(set_key(indi, "indi_pt"), on=["rid", "quarter"], how="left")
    c = c.merge(other_drug_keys(pull, rids_by_q), on=["rid", "quarter"], how="left")
    c[["pt_key", "indi_pt_key", "other_key"]] = c[["pt_key", "indi_pt_key", "other_key"]].fillna("")

    occ = c["occr_country"].astype(str).replace({"nan": ""})
    rep = c["reporter_country"].astype(str).replace({"nan": ""})
    country = occ.where(occ != "", rep)
    sex = c["sex"].astype(str).replace({"nan": ""})
    date = pd.Series(np.where(c["schema"] == "F", c["init_fda_dt"], c["fda_dt"]), index=c.index).astype("Int64")
    drug = pd.Series("", index=c.index)
    for i in ING:
        drug = drug + np.where(c[f"{i}_psss"], i + "|", "")
    key = pd.DataFrame({"age": c["age_years"].astype("float64").round(1), "sex": sex, "country": country,
                        "date": date, "drug": drug, "other": c["other_key"], "indi": c["indi_pt_key"],
                        "reac": c["pt_key"]})
    missing = key["age"].isna() | (key["sex"] == "") | (key["country"] == "")
    ks = key.iloc[np.argsort(c["caseid"].to_numpy(), kind="stable")]
    keep_a = pd.Series(True, index=key.index)
    keep_a.loc[ks.index[ks.duplicated(keep="first").to_numpy()]] = False
    kn = ks[~missing.loc[ks.index]]
    keep_b = pd.Series(True, index=key.index)
    keep_b.loc[kn.index[kn.duplicated(keep="first").to_numpy()]] = False
    size = key.assign(_n=1).groupby(list(key.columns), dropna=False)["_n"].transform("size")

    res = c[["rid", "quarter", "caseid"]].copy()
    res["keep_missing_as_value"] = keep_a.to_numpy()
    res["keep_conservative"] = keep_b.to_numpy()
    res["group_size"] = size.to_numpy()
    res.to_parquet(d / "fq_cases_dedup.parquet", index=False)

    summary = {"fq_cases": int(len(c)), "legacy_cases_using_fda_dt_as_date": int((c["schema"] == "L").sum()),
               "cases_missing_age_sex_or_country": int(missing.sum()),
               "removed_missing_as_value": int((~keep_a).sum()), "removed_conservative": int((~keep_b).sum()),
               "paper_table1": PAPER}
    for name, keep in (("missing_as_value", keep_a), ("conservative", keep_b)):
        summary[name] = summarize(c, keep)
    (d / "02b_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"\nFQ cases {len(c):,}; removed as duplicates: {int((~keep_a).sum()):,} (missing as value), "
          f"{int((~keep_b).sum()):,} (conservative)")
    cols_p = ("all", *ING, "combination")
    for variant in ("missing_as_value", "conservative"):
        print(f"\n{variant}:")
        for k in ("A_no_other_suspect", "B_no_other_drug_any_role"):
            for tag in ("before", "after"):
                r = summary[variant][f"{k}_{tag}"]
                print(f"  {k[:1]} {tag:6s}: " + ", ".join(f"{a}={r[a]:,}" for a in cols_p))
    print("  paper     : " + ", ".join(f"{a}={PAPER[a]:,}" for a in cols_p))
    print(f"\nWritten to {d}")


ROOT_DIR = s01.ROOT

if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in s01.PULLS:
        sys.exit("usage: 02b_deduplicate_fq_cases.py [validation|extension]")
    main(sys.argv[1])
