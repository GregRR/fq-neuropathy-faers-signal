"""
01_parse_and_link_core_files.py

Parses DEMO, DRUG and REAC from the FDA quarterly ASCII ZIPs, harmonizes the two
file layouts (legacy AERS, keyed on ISR; and FAERS, keyed on PRIMARYID), applies the
case-versioning dedup rule, and writes slim linked tables to data/interim/<pull>/.

Usage (from the Flox folder):
    python3 scripts/01_parse_and_link_core_files.py validation   # 2007q1 .. 2024q3
    python3 scripts/01_parse_and_link_core_files.py extension    # 2024q4 .. 2026q2
Requires: pandas, pyarrow   (pip install pandas pyarrow)

Layout facts used (FDA ASC_NTS.pdf section G, Readme.pdf, and header rows of the
actual files):
  * Legacy AERS (aers_ascii_*, through 2012q3): ISR, CASE, FOLL_SEQ, GNDR_COD, ...
    One CASE can have several ISRs; the latest ISR is the most current (Readme.pdf).
    DRUG/REAC keyed on ISR only. No INIT_FDA_DT, no ROLE-independent ingredient field.
  * FAERS (faers_ascii_*, 2012q4 onward): PRIMARYID, CASEID, CASEVERSION, INIT_FDA_DT;
    GNDR_COD renamed SEX from 2014q3; PROD_AI (active ingredient) added 2014q3.
    Each quarterly file holds only the latest version of a case available at
    extract time, so a case can appear in several quarters with rising CASEVERSION.
  * Deleted/DELETE*.txt (present in recent ZIPs) lists cases FDA removed.

DEDUP RULE (stated for METHODS_NOTES.md), applied within one pull:
  1. Drop any CASEID listed in a Deleted file of a quarter in the pull.
  2. Group by CASEID. FAERS-era rows outrank legacy rows. Within FAERS keep the highest
     CASEVERSION (ties: latest FDA_DT, then highest PRIMARYID). Within legacy keep the
     latest FDA_DT (ties: highest ISR).
  3. Keep exactly one report per CASEID. Its DRUG and REAC rows are the ones from the
     same quarter file as the retained report.
Each pull is deduplicated on its own, so the validation window sees only data FDA had
released through 2024q3, as the published paper did.

Outputs in data/interim/<pull>/:
  demo.parquet           one row per retained case
  drug/<quarter>.parquet rid, drug_seq, role_cod, drugname, prod_ai (upper-cased)
  reac/<quarter>.parquet rid, pt (lower-cased)
  01_summary.json        row counts, malformed-line counts, overlap and orphan counts
`rid` = 'L<ISR>' (legacy) or 'F<PRIMARYID>' (FAERS) so the two ID spaces cannot collide.
"""
import csv, json, re, sys, zipfile
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PULLS = {"validation": ((2007, 1), (2024, 3)), "extension": ((2024, 4), (2026, 2))}
AGE_YEARS = {"YR": 1.0, "DEC": 10.0, "MON": 1 / 12, "WK": 1 / 52, "DY": 1 / 365, "HR": 1 / 8760}
MEMBER = re.compile(r"^(demo|drug|reac|delete)\d{2}q[1-4]\.txt$", re.I)


def quarters(start, end):
    y, q = start
    while (y, q) <= end:
        yield y, q
        y, q = (y + 1, 1) if q == 4 else (y, q + 1)


def zip_path(y, q):
    return RAW / f"{'aers' if (y, q) <= (2012, 3) else 'faers'}_ascii_{y}q{q}.zip"


def members(zf):
    out = {}
    for n in zf.namelist():
        m = MEMBER.match(Path(n).name)
        if m:
            out[m.group(1).lower()] = n
    return out


def count_lines(zf, name):
    n = 0
    with zf.open(name) as fh:
        while chunk := fh.read(1 << 24):
            n += chunk.count(b"\n")
    return n


def read_member(zf, name, stats, key):
    with zf.open(name) as fh:
        df = pd.read_csv(fh, sep="$", dtype=str, quoting=csv.QUOTE_NONE, encoding="utf-8-sig",
                         encoding_errors="replace", index_col=False, on_bad_lines="skip",
                         keep_default_na=False)
    df.columns = [c.strip().lower() for c in df.columns]
    df = df.fillna("")
    skipped = count_lines(zf, name) - 1 - len(df)
    stats["lines_skipped"][key] = max(skipped, 0)
    return df


def tidy_demo(df, schema, quarter):
    out = pd.DataFrame(index=df.index)
    if schema == "L":
        out["primaryid"] = df["isr"]; out["caseid"] = df["case"]
        out["caseversion"] = "0"; out["init_fda_dt"] = ""
        out["sex"] = df["gndr_cod"]; out["occr_country"] = ""
    else:
        out["primaryid"] = df["primaryid"]; out["caseid"] = df["caseid"]
        out["caseversion"] = df["caseversion"]; out["init_fda_dt"] = df["init_fda_dt"]
        out["sex"] = df["sex"] if "sex" in df else df["gndr_cod"]
        out["occr_country"] = df.get("occr_country", "")
    out["fda_dt"] = df["fda_dt"]
    out["reporter_country"] = df.get("reporter_country", "")
    age = pd.to_numeric(df["age"], errors="coerce")
    out["age_years"] = (age * df["age_cod"].str.strip().str.upper().map(AGE_YEARS)).astype("float32")
    for c in ("primaryid", "caseid", "caseversion", "init_fda_dt", "fda_dt"):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("Int64")
    out["sex"] = out["sex"].str.strip().str.upper().astype("category")
    out["reporter_country"] = out["reporter_country"].str.strip().str.upper().astype("category")
    out["occr_country"] = out["occr_country"].str.strip().str.upper().astype("category")
    out["schema"] = schema
    out["quarter"] = quarter
    out["rid"] = schema + out["primaryid"].astype(str)
    return out


def read_deleted(zf, name, quarter, stats):
    raw = pd.read_csv(zf.open(name), sep="$", dtype=str, header=None, index_col=False,
                      quoting=csv.QUOTE_NONE, keep_default_na=False)
    first = raw.iloc[:, 0].str.strip()
    ids = first[first.str.fullmatch(r"\d+")]
    nonnum = len(first) - len(ids)
    print(f"  {quarter}: Deleted file, {len(ids)} numeric ids, {nonnum} non-numeric rows "
          f"(first rows: {list(first.head(3))})")
    stats["deleted_nonnumeric_rows"][quarter] = nonnum
    return set(ids.astype("int64"))


def main(pull):
    start, end = PULLS[pull]
    qs = list(quarters(start, end))
    outdir = ROOT / "data" / "interim" / pull
    (outdir / "drug").mkdir(parents=True, exist_ok=True)
    (outdir / "reac").mkdir(parents=True, exist_ok=True)
    stats = {"pull": pull, "quarters": len(qs), "lines_skipped": {}, "deleted_nonnumeric_rows": {}}

    # Pass 1: DEMO for every quarter, plus Deleted lists
    demos, deleted = [], set()
    for y, q in qs:
        label, zp = f"{y}q{q}", zip_path(y, q)
        if not zp.exists():
            sys.exit(f"Missing {zp}. Run 00_download_faers_quarters.py first.")
        with zipfile.ZipFile(zp) as zf:
            mem = members(zf)
            schema = "L" if (y, q) <= (2012, 3) else "F"
            demos.append(tidy_demo(read_member(zf, mem["demo"], stats, f"{label}:demo"), schema, label))
            if "delete" in mem:
                deleted |= read_deleted(zf, mem["delete"], label, stats)
        print(f"demo {label}: {len(demos[-1]):,} reports", flush=True)
    demo = pd.concat(demos, ignore_index=True)
    del demos
    stats["demo_rows_read"] = len(demo)

    # Dedup (rule in module docstring)
    n_del = int(demo["caseid"].isin(deleted).sum())
    demo = demo[~demo["caseid"].isin(deleted)]
    stats["reports_removed_by_deleted_lists"] = n_del
    both = set(demo.loc[demo.schema == "L", "caseid"]) & set(demo.loc[demo.schema == "F", "caseid"])
    stats["caseids_in_both_legacy_and_faers"] = len(both)
    demo = demo.assign(_rank=(demo.schema == "F").astype(int), _v=demo.caseversion.fillna(0),
                       _d=demo.fda_dt.fillna(0))
    demo = (demo.sort_values(["caseid", "_rank", "_v", "_d", "primaryid"])
                .drop_duplicates("caseid", keep="last")
                .drop(columns=["_rank", "_v", "_d"]).reset_index(drop=True))
    stats["cases_retained"] = len(demo)
    demo.to_parquet(outdir / "demo.parquet", index=False)
    print(f"retained {len(demo):,} cases (removed {n_del:,} reports via Deleted lists)")

    # Pass 2: DRUG and REAC, restricted to retained report of each case
    keep_by_q = {k: set(v) for k, v in demo.groupby("quarter")["rid"]}
    no_drug, no_reac = 0, 0
    for y, q in qs:
        label = f"{y}q{q}"
        schema = "L" if (y, q) <= (2012, 3) else "F"
        keep = keep_by_q.get(label, set())
        with zipfile.ZipFile(zip_path(y, q)) as zf:
            mem = members(zf)
            dr = read_member(zf, mem["drug"], stats, f"{label}:drug")
            if schema == "L":
                dr = dr.rename(columns={"isr": "primaryid"})
            dr["rid"] = schema + dr["primaryid"]
            dr = dr[dr["rid"].isin(keep)]
            if "prod_ai" not in dr: dr["prod_ai"] = ""
            dr = dr[["rid", "drug_seq", "role_cod", "drugname", "prod_ai"]].copy()
            dr["role_cod"] = dr["role_cod"].str.strip().str.upper()
            dr["drugname"] = dr["drugname"].str.strip().str.upper()
            dr["prod_ai"] = dr["prod_ai"].str.strip().str.upper()
            dr.to_parquet(outdir / "drug" / f"{label}.parquet", index=False)

            rc = read_member(zf, mem["reac"], stats, f"{label}:reac")
            if schema == "L":
                rc = rc.rename(columns={"isr": "primaryid"})
            rc["rid"] = schema + rc["primaryid"]
            rc = rc[rc["rid"].isin(keep)][["rid", "pt"]].copy()
            rc["pt"] = rc["pt"].str.strip().str.lower()
            rc.to_parquet(outdir / "reac" / f"{label}.parquet", index=False)
        no_drug += len(keep - set(dr["rid"])); no_reac += len(keep - set(rc["rid"]))
        print(f"link {label}: {len(keep):,} cases, {len(dr):,} drug rows, {len(rc):,} reaction rows", flush=True)

    stats["retained_cases_without_drug_rows"] = no_drug
    stats["retained_cases_without_reaction_rows"] = no_reac
    (outdir / "01_summary.json").write_text(json.dumps(stats, indent=2, default=int))
    print(f"\nDone. Summary: {outdir / '01_summary.json'}")


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in PULLS:
        sys.exit("usage: 01_parse_and_link_core_files.py [validation|extension]")
    main(sys.argv[1])
