"""
08_comparator_controls.py

Control drugs and comparator antibiotics for the fluoroquinolone (FQ) peripheral neuropathy (PN) signal, run through
the same case rules, duplicate removal, event definition and backgrounds as the FQ groups.

Guo et al. 2026 (Int J Antimicrob Agents 67:107993) used azithromycin as a negative control (no signal expected) and
metronidazole as a positive control (neuropathy is an established adverse effect). Their drug definitions are in a
supplementary table we do not have, so the definitions below are our own.

The comparator antibiotics (amoxicillin with or without clavulanate, cephalexin, sulfamethoxazole-trimethoprim, and
azithromycin, which is also the negative control) are the comparison drugs of script 09. Rule A and B counts for them
are built here and used there.

Usage (from the repository root, after scripts 01, 02, 02b, 03, 04 and 05):
    python3 scripts/08_comparator_controls.py validation
    python3 scripts/08_comparator_controls.py extension
    python3 scripts/08_comparator_controls.py summary      # after both pulls

Drug identification, as for the FQs in script 02: DRUGNAME and PROD_AI are searched for generic-name stems and
U.S. brand names; names the FDA's PROD_AI links to exactly one group, but the patterns miss, are added from a
lookup built over both pulls (needed for legacy files, which have no PROD_AI).
  azithromycin   azithromycin, Zithromax (including Z-Pak and TRI-PAK), Zmax, Azasite
  metronidazole  metronidazole, Flagyl, Likmez, Metro I.V., Protostat, Metryl, MetroCream, MetroGel, MetroLotion,
                 Noritate, Rosadan, Nuvessa, Vandazole
  amoxicillin    amoxicillin, clavulanate or clavulanic acid, Amoxil, Trimox, Moxilin, Moxatag, Augmentin
                 (including XR and ES-600), Amoclan. Clavulanate rows are counted in this group so that a product
                 reported as two rows is not treated as having another drug.
  cephalexin     cephalexin, cefalexin, Keflex, Biocef, Zartan, Panixine, Daxbia
  sulfamethoxazole_trimethoprim  sulfamethoxazole, trimethoprim, co-trimoxazole, Bactrim, Septra, Sulfatrim, Bethaprim,
                 Cotrim, Uroplus, SMZ-TMP. Trimethoprim alone is included.
Brand names: drugs.com pages for each drug (the metronidazole and amoxicillin pages hide some brands, so those lists
may be incomplete). The multi-agent Helicobacter pylori products
Pylera (bismuth subcitrate potassium, metronidazole, tetracycline), Prevpac (amoxicillin, clarithromycin,
lansoprazole), Talicia (amoxicillin, omeprazole, rifabutin), Omeclamox-Pak (amoxicillin, clarithromycin, omeprazole)
and the Voquezna packs (amoxicillin with vonoprazan, with or without clarithromycin) are excluded from every group;
ingredients are from the drugs.com H. pylori eradication agents class page.
Topical products (ophthalmic azithromycin; topical and vaginal metronidazole) are flagged by route or name, and the
scenario 'excl_topical_only' drops cases whose suspect rows are all topical.

Case rules, as for the FQs (script 02), with 'other drug' meaning any drug outside the group:
  A  group drug as primary or secondary suspect, no other drug as primary or secondary suspect
  B  group drug as primary or secondary suspect, no other drug in any role
Duplicate removal follows script 02b (age, sex, country, first FDA date, suspect drug, other drugs, indications,
reactions; lowest CASEID kept), both variants. Event: any PN24 term (script 03).

Backgrounds are those of script 04 (backgrounds.csv). Rule B cases contain no other drug, so they lie inside every
background (a case with the group drug in a second, non-suspect row is the rare exception). Rule A is reported on
the 'all' and 'paper_implied' backgrounds only, the two that contain every report.

Output in data/interim/<pull>/:
  comparator_cases.parquet   one row per case and group: rid, quarter, caseid, group, rule flags, topical-only,
                             pn24, duplicate-removal flags
  comparator_signals.csv     counts, ROR and IC for each group, rule, duplicate variant, topical scenario, background
  comparator_name_audit.csv  every DRUGNAME and PROD_AI pair matched to each group, with row counts (review for
                             combination products that should not be in a group)
  comparator_name_map.csv    names added from the PROD_AI lookup
Output of 'summary':
  results/tables/comparator_controls.csv   rule B, variant 1, all cases, PN24: controls beside the FQ groups, both pulls
  results/figures/comparator_controls.png, .pdf   IC with 95% interval on two backgrounds, validation pull
"""
import importlib.util
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s02b = load("s02b", "02b_deduplicate_fq_cases.py")
s03b = load("s03b", "03b_flag_topical_cases.py")
s03 = load("s03", "03_filter_neuropathy_reactions.py")
s05 = load("s05", "05_compute_ror_and_bcpnn.py")
ROOT = s02b.ROOT_DIR

GROUPS = {
    "azithromycin": r"\bAZITHROMYCIN|\bZITHROMAX\b|\bZ-?PAK\b|\bZMAX\b|\bAZASITE\b",
    "metronidazole": r"\bMETRONIDAZOL|\bFLAGYL\b|\bLIKMEZ\b|\bMETRO I\.?V\b|\bPROTOSTAT\b|\bMETRYL\b|\bMETROCREAM\b|"
                     r"\bMETROGEL\b|\bMETROLOTION\b|\bNORITATE\b|\bROSADAN\b|\bNUVESSA\b|\bVANDAZOLE\b",
    "amoxicillin": r"\bAMOXICILLIN|CLAVULAN|\bAMOXIL\b|\bTRIMOX\b|\bMOXILIN\b|\bMOXATAG\b|\bAUGMENTIN\b|\bAMOCLAN\b",
    "cephalexin": r"\bCEPHALEXIN|\bCEFALEXIN|\bKEFLEX\b|\bBIOCEF\b|\bZARTAN\b|\bPANIXINE\b|\bDAXBIA\b",
    "sulfamethoxazole_trimethoprim": r"\bSULFAMETHOXAZOL|\bTRIMETHOPRIM|\bCO-?TRIMOX|\bBACTRIM|\bSEPTRA\b|\bSULFATRIM|"
                                     r"\bBETHAPRIM|\bCOTRIM\b|\bUROPLUS|\bSMZ[- ]?TMP",
}
RX = {g: re.compile(v) for g, v in GROUPS.items()}
COMBO = re.compile(r"\bPYLERA\b|\bPREVPAC\b|\bTALICIA\b|\bOMECLAMOX|\bVOQUEZNA\b")
MIN_ROWS = 3
TOPICAL_ROUTE = set(s03b.TOPICAL_ROUTE) | {"VAGINAL"}
TOPICAL_NAME = re.compile(r"OPHTH|\bOTIC\b|\bEYES?\b|\bEARS?\b|DROP|OINTMENT|CREAM|\bGEL\b|LOTION|VAGINAL|\bTOPICAL\b|"
                          r"\bAZASITE\b|\bMETROCREAM\b|\bMETROGEL\b|\bMETROLOTION\b|\bNORITATE\b|\bROSADAN\b|"
                          r"\bNUVESSA\b|\bVANDAZOLE\b")
RULES = {"A": "A_no_other_suspect", "B": "B_no_other_drug_any_role"}
CONTROLS = ("azithromycin", "metronidazole")
S24 = set(s03.PN24)


def classify(uniq, name_map=None):
    """bool column per group for the distinct (drugname, prod_ai) pairs in uniq."""
    out = pd.DataFrame(index=uniq.index)
    combo = uniq["drugname"].str.contains(COMBO)
    for g, rx in RX.items():
        f = uniq["drugname"].str.contains(rx) | uniq["prod_ai"].str.contains(rx)
        if name_map:
            f |= uniq["drugname"].isin(name_map.get(g, ()))
        out[g] = f & ~combo
    return out


def build_name_map(parts):
    """Names the FDA's PROD_AI links to exactly one group but the patterns do not match."""
    frames = []
    for p in parts:
        df = pd.read_parquet(p, columns=["drugname", "prod_ai"]).fillna("")
        df = df[(df["prod_ai"] != "") & (df["drugname"] != "")]
        g = df.groupby(["drugname", "prod_ai"]).size().rename("rows").reset_index()
        for grp, rx in RX.items():
            sel = g[g["prod_ai"].str.contains(rx) & ~g["drugname"].str.contains(rx) & ~g["drugname"].str.contains(COMBO)]
            if len(sel):
                frames.append(sel.assign(group=grp)[["drugname", "prod_ai", "group", "rows"]])
    if not frames:
        return {}, pd.DataFrame(columns=["drugname", "prod_ai", "group", "rows"])
    tab = pd.concat(frames).groupby(["drugname", "prod_ai", "group"], as_index=False)["rows"].sum()
    n_grp = tab.groupby("drugname")["group"].transform("nunique")
    keep = tab[(tab["rows"] >= MIN_ROWS) & (n_grp == 1)].sort_values("rows", ascending=False)
    return {grp: set(g["drugname"]) for grp, g in keep.groupby("group")}, keep


def process_partition(path, name_map, name_counts):
    df = pd.read_parquet(path)
    for c in ("drugname", "prod_ai", "route"):
        df[c] = df[c].fillna("")
    uniq = df[["drugname", "prod_ai"]].drop_duplicates().reset_index(drop=True)
    flags = classify(uniq, name_map)
    df = df.merge(pd.concat([uniq, flags], axis=1), on=["drugname", "prod_ai"], how="left")
    ps = df["role_cod"].isin(["PS", "SS"])
    topical = df["route"].isin(TOPICAL_ROUTE) | df["drugname"].str.contains(TOPICAL_NAME)
    cols = []
    for g in GROUPS:
        df[f"{g}_any"] = df[g]
        df[f"{g}_psss"] = df[g] & ps
        df[f"{g}_nt"] = df[g] & ps & ~topical
        df[f"{g}_oth_any"] = ~df[g]
        df[f"{g}_oth_ps"] = ~df[g] & ps
        cols += [f"{g}_any", f"{g}_psss", f"{g}_nt", f"{g}_oth_any", f"{g}_oth_ps"]
        name_counts.append(df[df[g]].groupby(["drugname", "prod_ai"]).size().rename("rows").reset_index().assign(group=g))
    cases = df.groupby("rid", sort=False)[cols].max()
    cases = cases[cases[[f"{g}_any" for g in GROUPS]].any(axis=1)].reset_index()
    cases.insert(1, "quarter", path.stem)
    # drugs outside each group on the cases with a suspect group drug, for the duplicate key
    ident = df["prod_ai"].where(df["prod_ai"] != "", df["drugname"])
    others = []
    for g in GROUPS:
        cand = set(cases.loc[cases[f"{g}_psss"], "rid"])
        m = df["rid"].isin(cand) & ~df[g] & (ident != "")
        others.append(pd.DataFrame({"rid": df.loc[m, "rid"], "ident": ident[m]}).assign(group=g, quarter=path.stem))
    return cases, pd.concat(others, ignore_index=True)


def dedup_flags(key, caseid):
    """Both keep flags, as in script 02b (lowest CASEID kept within a group of identical keys)."""
    missing = key["age"].isna() | (key["sex"] == "") | (key["country"] == "")
    ks = key.iloc[np.argsort(caseid.to_numpy(), kind="stable")]
    keep_a = pd.Series(True, index=key.index)
    keep_a.loc[ks.index[ks.duplicated(keep="first").to_numpy()]] = False
    kn = ks[~missing.loc[ks.index]]
    keep_b = pd.Series(True, index=key.index)
    keep_b.loc[kn.index[kn.duplicated(keep="first").to_numpy()]] = False
    return keep_a, keep_b


def build_cases(pull, cases, others):
    d = ROOT / "data" / "interim" / pull
    demo = pd.read_parquet(d / "demo.parquet")
    dcols = ["rid", "caseid", "quarter", "schema", "age_years", "sex", "reporter_country", "occr_country",
             "init_fda_dt", "fda_dt"]
    rids_by_q = {}
    for g in GROUPS:
        for q, s in cases[cases[f"{g}_psss"]].groupby("quarter")["rid"]:
            rids_by_q.setdefault(q, set()).update(s)
    reac = []
    for q, keep in rids_by_q.items():
        r = pd.read_parquet(d / "reac" / f"{q}.parquet", columns=["rid", "pt"])
        reac.append(r[r["rid"].isin(keep)].assign(quarter=q))
    reac = pd.concat(reac, ignore_index=True)
    reac["pn"] = s03.normalise(reac["pt"])
    pn24 = reac[reac["pn"].isin(S24)][["rid", "quarter"]].drop_duplicates().assign(pn24=True)
    indi = s02b.read_indications(pull, rids_by_q)
    pt_key = s02b.set_key(reac, "pt")
    indi_key = s02b.set_key(indi, "indi_pt")

    out = []
    for g in GROUPS:
        sub = cases[cases[f"{g}_psss"]].merge(demo[dcols], on=["rid", "quarter"], how="left", validate="one_to_one")
        o = (others[others["group"] == g].groupby(["rid", "quarter"])["ident"]
             .agg(lambda s: "|".join(sorted(set(s)))).rename("other_key").reset_index())
        sub = (sub.merge(pt_key, on=["rid", "quarter"], how="left").merge(indi_key, on=["rid", "quarter"], how="left")
               .merge(o, on=["rid", "quarter"], how="left").merge(pn24, on=["rid", "quarter"], how="left"))
        sub[["pt_key", "indi_pt_key", "other_key"]] = sub[["pt_key", "indi_pt_key", "other_key"]].fillna("")
        sub["pn24"] = sub["pn24"].fillna(False).astype(bool)
        occ = sub["occr_country"].astype(str).replace({"nan": ""})
        rep = sub["reporter_country"].astype(str).replace({"nan": ""})
        key = pd.DataFrame({
            "age": sub["age_years"].astype("float64").round(1),
            "sex": sub["sex"].astype(str).replace({"nan": ""}),
            "country": occ.where(occ != "", rep),
            "date": pd.Series(np.where(sub["schema"] == "F", sub["init_fda_dt"], sub["fda_dt"]), index=sub.index).astype("Int64"),
            "drug": g, "other": sub["other_key"], "indi": sub["indi_pt_key"], "reac": sub["pt_key"]})
        ka, kb = dedup_flags(key, sub["caseid"])
        out.append(pd.DataFrame({
            "rid": sub["rid"], "quarter": sub["quarter"], "caseid": sub["caseid"], "group": g,
            "rule_A": ~sub[f"{g}_oth_ps"], "rule_B": ~sub[f"{g}_oth_any"],
            "topical_only": ~sub[f"{g}_nt"], "pn24": sub["pn24"],
            "keep_missing_as_value": ka.to_numpy(), "keep_conservative": kb.to_numpy()}))
    return pd.concat(out, ignore_index=True)


def signals(cases, bg):
    rng = np.random.default_rng(s05.SEED)
    rows = []
    for g in GROUPS:
        cs = cases[cases["group"] == g]
        for rule in RULES:
            for var, keep in (("missing_as_value", cs["keep_missing_as_value"]), ("conservative", cs["keep_conservative"])):
                for scen, sm in (("all_cases", pd.Series(True, index=cs.index)), ("excl_topical_only", ~cs["topical_only"])):
                    k = cs[f"rule_{rule}"] & keep & sm
                    nd, o = int(k.sum()), int((k & cs["pn24"]).sum())
                    for b, r in bg.iterrows():
                        if rule == "A" and b not in ("all", "paper_implied"):
                            continue
                        ne, n = int(r["n_pn24"]), int(r["N"])
                        row = {"group": g, "rule": RULES[rule], "dedup": var, "topical": scen, "event": "pn24",
                               "background": b, "O": o, "n_drug": nd, "n_event": ne, "N": n,
                               "a": o, "b": nd - o, "c": ne - o, "d": n - nd - (ne - o)}
                        est = dict.fromkeys(["ror", "ror_lo", "ror_hi", "p_fisher", "log10_p_fisher", "ic", "ic025", "ic975"], np.nan)
                        if o > 0 and nd > 0 and min(row["a"], row["b"], row["c"], row["d"]) >= 0:
                            e = s05.ror(row["a"], row["b"], row["c"], row["d"])
                            p = s05.fisher_p(row["a"], row["b"], row["c"], row["d"])
                            ic = s05.ic_posterior(o, nd, ne, n, rng)
                            est = dict(zip(est, (*e, *p, *ic)))
                        rows.append({**row, **est})
    return pd.DataFrame(rows)


def main(pull):
    d = ROOT / "data" / "interim" / pull
    for f in ("demo.parquet", "backgrounds.csv"):
        if not (d / f).exists():
            sys.exit(f"Missing {d / f}. Run scripts 01 and 04 first.")
    parts = sorted((d / "drug").glob("*.parquet"))
    lookup = sorted(q for pl in ("validation", "extension") for q in (ROOT / "data" / "interim" / pl / "drug").glob("*.parquet"))
    name_map, name_tab = build_name_map(lookup)
    name_tab.to_csv(d / "comparator_name_map.csv", index=False)
    name_counts, cases, others = [], [], []
    for p in parts:
        c, o = process_partition(p, name_map, name_counts)
        cases.append(c)
        others.append(o)
        print(f"{p.stem}: " + ", ".join(f"{g} {int(c[f'{g}_any'].sum()):,}" for g in GROUPS), flush=True)
    cases, others = pd.concat(cases, ignore_index=True), pd.concat(others, ignore_index=True)
    pd.concat(name_counts).groupby(["group", "drugname", "prod_ai"], as_index=False)["rows"].sum().sort_values(
        ["group", "rows"], ascending=False).to_csv(d / "comparator_name_audit.csv", index=False)

    long = build_cases(pull, cases, others)
    long.to_parquet(d / "comparator_cases.parquet", index=False)
    bg = pd.read_csv(d / "backgrounds.csv").set_index("background").rename(columns={"n_pn24": "n_pn24"})
    res = signals(long, bg)
    res.to_csv(d / "comparator_signals.csv", index=False)

    sel = res[(res["rule"] == RULES["B"]) & (res["dedup"] == "missing_as_value") & (res["topical"] == "all_cases")]
    print("\nRule B, duplicate variant 1, all cases, PN24 event:")
    for b, g in sel.groupby("background"):
        print(f"\n  background {b}")
        for r in g.itertuples():
            if pd.isna(r.ic):
                print(f"    {r.group:14s} {r.O:>6,}/{r.n_drug:<7,}  no estimate")
            else:
                print(f"    {r.group:14s} {r.O:>6,}/{r.n_drug:<7,} IC {r.ic:5.2f} ({r.ic025:5.2f}, {r.ic975:5.2f})  "
                      f"ROR {r.ror:6.2f} ({r.ror_lo:6.2f}, {r.ror_hi:6.2f})  p {s05.fmt_p(r.log10_p_fisher)}")
    print(f"\n{len(res):,} rows written to {d / 'comparator_signals.csv'}")


FQ_SHOW = ("ciprofloxacin", "levofloxacin", "moxifloxacin", "all")


def summary():
    frames = []
    for pull in ("validation", "extension"):
        d = ROOT / "data" / "interim" / pull
        fq, cc = d / "signals.csv", d / "comparator_signals.csv"
        if not (fq.exists() and cc.exists()):
            sys.exit(f"Missing {fq} or {cc}. Run scripts 05 and 08 on the {pull} pull first.")
        f = pd.read_csv(fq)
        f = f[(f["event"] == "pn24") & (f["topical"] == "all_cases") & (f["group"].isin(FQ_SHOW))].assign(kind="fluoroquinolone")
        c = pd.read_csv(cc).assign(kind="control")
        t = pd.concat([f, c], ignore_index=True)
        t = t[(t["rule"] == RULES["B"]) & (t["dedup"] == "missing_as_value") & (t["topical"] == "all_cases")]
        frames.append(t.assign(pull=pull))
    t = pd.concat(frames, ignore_index=True)
    keep = ["pull", "kind", "group", "background", "O", "n_drug", "n_event", "N", "ror", "ror_lo", "ror_hi",
            "log10_p_fisher", "ic", "ic025", "ic975"]
    tdir = ROOT / "results" / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    t[keep].to_csv(tdir / "comparator_controls.csv", index=False, float_format="%.4f")

    order = [*CONTROLS, *FQ_SHOW]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    y0 = np.arange(len(order))[::-1].astype(float)
    for ax, bg in zip(axes, ("all", "no_concomitant")):
        s = t[(t["pull"] == "validation") & (t["background"] == bg)].set_index("group").reindex(order)
        ok = s["ic"].notna().to_numpy()
        col = ["tab:green" if g in CONTROLS else "tab:blue" for g in order]
        for i in np.flatnonzero(ok):
            ax.errorbar(s["ic"].iloc[i], y0[i], xerr=[[s["ic"].iloc[i] - s["ic025"].iloc[i]], [s["ic975"].iloc[i] - s["ic"].iloc[i]]],
                        fmt="o", color=col[i], ms=5, lw=1.2, capsize=2)
        for i in np.flatnonzero(~ok):
            ax.text(0.05, y0[i], "no estimate", va="center", fontsize=8)
        ax.axvline(0, color="grey", lw=0.8, ls=":")
        ax.axvline(2, color="grey", lw=0.8, ls="--")
        ax.set_title(f"Background '{bg}'", fontsize=10)
        ax.set_xlabel("Information Component (95% interval)")
    axes[0].set_yticks(y0)
    axes[0].set_yticklabels(order)
    fig.suptitle("Control drugs (green) beside fluoroquinolones (blue): validation pull, rule B", fontsize=10)
    fig.tight_layout()
    fdir = ROOT / "results" / "figures"
    fdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fdir / "comparator_controls.png", dpi=300)
    fig.savefig(fdir / "comparator_controls.pdf")
    plt.close(fig)

    pd.set_option("display.width", 200)
    for bg in ("all", "no_concomitant"):
        print(f"\nBackground '{bg}', rule B, variant 1, all cases")
        s = t[t["background"] == bg].pivot(index="group", columns="pull", values=["O", "n_drug", "ic", "ic025"]).reindex(order)
        s.columns = [f"{m}_{p[:3]}" for m, p in s.columns]
        print(s.round(2).to_string())
    print(f"\nTable in {tdir}, figure in {fdir}")


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) == 2 else ""
    if arg in s03.PULLS:
        main(arg)
    elif arg == "summary":
        summary()
    else:
        sys.exit("usage: 08_comparator_controls.py [validation|extension|summary]")
