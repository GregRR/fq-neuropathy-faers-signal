"""
09_fq_vs_comparators.py

Fluoroquinolone (FQ) peripheral neuropathy (PN) reporting against comparator antibiotics instead of all of FAERS.
Each FQ group is compared with each comparator antibiotic on its own and with the four comparators pooled, so the
dependence of the result on the comparator is visible.

Usage (from the repository root, after scripts 04 and 08 on both pulls):
    python3 scripts/09_fq_vs_comparators.py

Comparators (script 08): amoxicillin with or without clavulanate, cephalexin, sulfamethoxazole-trimethoprim and
azithromycin. Counts are those of the case-by-case rules of scripts 02 to 04 and 08: rule A or B, duplicate variant 1
(missing_as_value), all cases or excluding topical-only cases, PN24 event. The background of each comparison is the
FQ group plus the comparator reports, so for FQ group g and comparator c
    a = PN reports with g           b = other reports with g
    c = PN reports with c           d = other reports with c
ROR with a Woolf interval and Fisher exact p-value, and the Information Component (Noren et al. 2006) with
n.. = n_g + n_c, n_event = O_g + O_c, as in script 05. The two pulls are disjoint, so the pooled rows add their counts.

The comparators are chosen drugs, not all other antibiotics. Nitrofurantoin, metronidazole and linezolid have
reported neuropathy associations and are not used as comparators.

Output:
  results/tables/fq_vs_comparators.csv   pull (validation, extension, pooled), rule, topical scenario, FQ group,
                                         comparator, counts, ROR, IC
  results/figures/fq_vs_comparators.png, .pdf   ROR with 95% interval on a log scale, validation pull, rule B
"""
import importlib.util
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
spec = importlib.util.spec_from_file_location("s05", HERE / "05_compute_ror_and_bcpnn.py")
s05 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s05)

RULES = {"A": "A_no_other_suspect", "B": "B_no_other_drug_any_role"}
TOPICAL = ("all_cases", "excl_topical_only")
FQ_GROUPS = ("all", "ciprofloxacin", "levofloxacin", "moxifloxacin", "combination", "ofloxacin", "gemifloxacin",
             "delafloxacin")
COMPARATORS = ("amoxicillin", "cephalexin", "sulfamethoxazole_trimethoprim", "azithromycin")
POOLED_COMPARATOR = "all_four"


def counts(path, groups):
    """O and n_drug per (rule, topical, group) for duplicate variant 1, PN24 event, background 'all'."""
    t = pd.read_csv(path)
    t = t[(t["dedup"] == "missing_as_value") & (t["event"] == "pn24") & (t["background"] == "all") &
          (t["group"].isin(groups))]
    return t.set_index(["rule", "topical", "group"])[["O", "n_drug"]]


def load_pull(pull):
    d = ROOT / "data" / "interim" / pull
    f, c = d / "contingency.csv", d / "comparator_signals.csv"
    for p in (f, c):
        if not p.exists():
            sys.exit(f"Missing {p}. Run scripts 04 and 08 on the {pull} pull first.")
    return counts(f, FQ_GROUPS), counts(c, COMPARATORS)


def compare(pull, fq, cmp_, rng):
    rows = []
    for rule, rname in RULES.items():
        for scen in TOPICAL:
            comp = {c: cmp_.loc[(rname, scen, c)] if (rname, scen, c) in cmp_.index else None for c in COMPARATORS}
            comp[POOLED_COMPARATOR] = pd.Series({"O": sum(v["O"] for v in comp.values() if v is not None),
                                                 "n_drug": sum(v["n_drug"] for v in comp.values() if v is not None)})
            for g in FQ_GROUPS:
                if (rname, scen, g) not in fq.index:
                    continue
                of, nf = int(fq.loc[(rname, scen, g), "O"]), int(fq.loc[(rname, scen, g), "n_drug"])
                for c, v in comp.items():
                    if v is None:
                        continue
                    oc, nc = int(v["O"]), int(v["n_drug"])
                    row = {"pull": pull, "rule": rule, "topical": scen, "fq_group": g, "comparator": c,
                           "O_fq": of, "n_fq": nf, "O_cmp": oc, "n_cmp": nc,
                           "share_fq": of / nf if nf else np.nan, "share_cmp": oc / nc if nc else np.nan}
                    est = dict.fromkeys(["ror", "ror_lo", "ror_hi", "log10_p_fisher", "ic", "ic025", "ic975"], np.nan)
                    a, b, c_, d = of, nf - of, oc, nc - oc
                    if of > 0 and nf > 0 and nc > 0:
                        e = s05.ror(a, b, c_, d)
                        _, lg = s05.fisher_p(a, b, c_, d)
                        ic = s05.ic_posterior(of, nf, of + oc, nf + nc, rng)
                        est = dict(zip(est, (*e, lg, *ic)))
                    rows.append({**row, **est})
    return pd.DataFrame(rows)


def pooled_counts(a, b):
    """Add the counts of two disjoint pulls."""
    idx = a.index.union(b.index)
    return a.reindex(idx).fillna(0).add(b.reindex(idx).fillna(0))


def forest(res, out):
    d = res[(res["pull"] == "validation") & (res["rule"] == "B") & (res["topical"] == "all_cases")]
    series = (*COMPARATORS, POOLED_COMPARATOR)
    labels = {"amoxicillin": "Amoxicillin (+/- clavulanate)", "cephalexin": "Cephalexin",
              "sulfamethoxazole_trimethoprim": "Sulfamethoxazole-trimethoprim", "azithromycin": "Azithromycin",
              POOLED_COMPARATOR: "Four pooled"}
    colors = ("tab:blue", "tab:orange", "tab:green", "tab:red", "black")
    marks = ("o", "s", "^", "D", "P")
    y0 = np.arange(len(FQ_GROUPS))[::-1].astype(float)
    offs = np.linspace(0.3, -0.3, len(series))
    fig, ax = plt.subplots(figsize=(8, 7))
    for c, col, mk, o in zip(series, colors, marks, offs):
        s = d[d["comparator"] == c].set_index("fq_group").reindex(list(FQ_GROUPS))
        ok = s["ror"].notna().to_numpy()
        ax.errorbar(s["ror"][ok], y0[ok] + o, xerr=[(s["ror"] - s["ror_lo"])[ok], (s["ror_hi"] - s["ror"])[ok]],
                    fmt=mk, color=col, ms=4, lw=1.1, capsize=2, label=labels[c])
    ax.axvline(1, color="grey", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_yticks(y0)
    ax.set_yticklabels(FQ_GROUPS)
    ax.set_xlabel("Reporting odds ratio against the comparator (95% interval, log scale)")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.set_title("Fluoroquinolone neuropathy reports against comparator antibiotics: validation pull, rule B", fontsize=10)
    fig.tight_layout()
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "fq_vs_comparators.png", dpi=300)
    fig.savefig(out / "fq_vs_comparators.pdf")
    plt.close(fig)


def main():
    fv, cv = load_pull("validation")
    fe, ce = load_pull("extension")
    rng = np.random.default_rng(s05.SEED)
    res = pd.concat([compare("validation", fv, cv, rng), compare("extension", fe, ce, rng),
                     compare("pooled", pooled_counts(fv, fe), pooled_counts(cv, ce), rng)], ignore_index=True)
    tdir = ROOT / "results" / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    res.to_csv(tdir / "fq_vs_comparators.csv", index=False, float_format="%.4f")
    forest(res, ROOT / "results" / "figures")

    pd.set_option("display.width", 220)
    for pull in ("validation", "extension"):
        s = res[(res["pull"] == pull) & (res["rule"] == "B") & (res["topical"] == "all_cases")]
        cs = s[s["fq_group"] == "all"].set_index("comparator")
        print(f"\n{pull} pull, rule B: neuropathy share (PN24 reports / reports) in each comparator")
        print("  " + ", ".join(f"{c} {cs.loc[c, 'O_cmp']:,.0f}/{cs.loc[c, 'n_cmp']:,.0f} ({100 * cs.loc[c, 'share_cmp']:.1f}%)"
                               for c in (*COMPARATORS, POOLED_COMPARATOR) if c in cs.index))
        print(f"\n{pull} pull, rule B: ROR (95% interval) of each FQ group against each comparator")
        w = s.assign(txt=s.apply(lambda r: "no estimate" if pd.isna(r["ror"]) else
                                 f"{r['ror']:.2f} ({r['ror_lo']:.2f}, {r['ror_hi']:.2f})", axis=1))
        w = w.pivot(index="fq_group", columns="comparator", values="txt").reindex(list(FQ_GROUPS))
        w = w[[c for c in (*COMPARATORS, POOLED_COMPARATOR) if c in w.columns]]
        w.columns = [c[:14] for c in w.columns]
        print(w.to_string())
    print(f"\nTable in {tdir}, figure in {ROOT / 'results' / 'figures'}")


if __name__ == "__main__":
    main()
