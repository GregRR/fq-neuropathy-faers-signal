"""
06_validate_against_table1.py

Compares the validation pull (script 05) with Table 1 of Shamim et al. 2026 (Clin Transl Sci 19:e70541).

Usage (from the repository root, after script 05 on the validation pull):
    python3 scripts/06_validate_against_table1.py

Selection from signals.csv: duplicate variant 1 (missing_as_value), all cases, PN24 event, for both sole-suspect
rules (A and B). Table 1 lists these groups: all fluoroquinolones, ciprofloxacin, levofloxacin, moxifloxacin,
combination, ofloxacin, gemifloxacin, delafloxacin.

Decomposition of the difference between our Information Component (IC) and the published one, per group:
  published IC                                                  (Table 1)
  IC from the published counts, paper-implied background        (recomputed; checks the estimator)
  IC from our counts, paper-implied background                  (the case-count effect is this minus the line above)
  IC from our counts, our background                            (the comparator effect is this minus the line above)
The paper-implied background is N = 17,610,899 reports with an event fraction of 1.826%, the values that
reproduce Table 1 (script 05, selftest).

Output:
  results/tables/validation_table.csv          our estimates next to Table 1, rule x background
  results/tables/validation_decomposition.csv  the decomposition above, rule x background
  results/figures/validation_forest.png, .pdf  IC with 95% interval: Table 1, our counts on our 'all' background,
                                               our counts on the paper-implied background (rule B)
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

N_PAPER = 17_610_899
N_EVENT_PAPER = round(0.01826 * N_PAPER)

# Table 1 of Shamim et al. 2026: reports with the drug and PN, reports with the drug, IC, lower and upper limit
PAPER = {"all": (5588, 37778, 3.02, 2.98, 3.05),
         "ciprofloxacin": (2291, 14434, 3.12, 3.06, 3.17),
         "levofloxacin": (1942, 13598, 2.97, 2.91, 3.02),
         "moxifloxacin": (642, 7439, 2.24, 2.13, 2.34),
         "combination": (652, 1505, 4.54, 4.46, 4.62),
         "ofloxacin": (46, 588, 2.05, 1.62, 2.42),
         "gemifloxacin": (11, 142, 1.89, 0.95, 2.59),
         "delafloxacin": (4, 72, 1.31, -0.39, 2.36)}
RULES = {"A": "A_no_other_suspect", "B": "B_no_other_drug_any_role"}
BACKGROUNDS = ("all", "excl_fq_any", "no_concomitant", "paper_implied")


def pick(t, rule, bg):
    s = t[(t["rule"] == RULES[rule]) & (t["dedup"] == "missing_as_value") & (t["topical"] == "all_cases") &
          (t["event"] == "pn24") & (t["background"] == bg)]
    return s.set_index("group").reindex(list(PAPER))


def build(t):
    tab, dec = [], []
    recalc = {g: s05.ic_posterior(o, nd, N_EVENT_PAPER, N_PAPER)[0] for g, (o, nd, *_) in PAPER.items()}
    for rule in RULES:
        pi = pick(t, rule, "paper_implied")
        for bg in BACKGROUNDS:
            s = pick(t, rule, bg)
            for g, (po, pnd, pic, plo, phi) in PAPER.items():
                r = s.loc[g]
                tab.append({"rule": rule, "background": bg, "group": g,
                            "paper_O": po, "paper_n_drug": pnd, "paper_ic": pic, "paper_ic_lo": plo, "paper_ic_hi": phi,
                            "O": r["O"], "n_drug": r["n_drug"], "ic": r["ic"], "ic025": r["ic025"], "ic975": r["ic975"],
                            "d_O": r["O"] - po, "d_n_drug": r["n_drug"] - pnd, "d_ic": r["ic"] - pic,
                            "paper_ic025_gt0": plo > 0, "ic025_gt0": bool(r["ic025"] > 0) if pd.notna(r["ic025"]) else np.nan})
                if bg in ("paper_implied",):
                    continue
                dec.append({"rule": rule, "background": bg, "group": g,
                            "paper_ic": pic, "paper_ic_recomputed": recalc[g],
                            "ours_ic_paper_background": pi.loc[g, "ic"], "ours_ic": r["ic"],
                            "case_count_effect": pi.loc[g, "ic"] - recalc[g],
                            "comparator_effect": r["ic"] - pi.loc[g, "ic"],
                            "rounding": recalc[g] - pic,
                            "total_gap": r["ic"] - pic})
    return pd.DataFrame(tab), pd.DataFrame(dec)


def forest(tab, out):
    d = tab[tab["rule"] == "B"]
    series = (("Shamim et al. 2026, Table 1", None, "black", "o"),
              ("Our counts, background 'all'", "all", "tab:blue", "s"),
              ("Our counts, paper-implied background", "paper_implied", "tab:orange", "^"))
    groups = list(PAPER)
    y0 = np.arange(len(groups))[::-1].astype(float)
    off = (0.22, 0.0, -0.22)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for (label, bg, col, mk), o in zip(series, off):
        if bg is None:
            s = d[d["background"] == "all"].set_index("group").reindex(groups)
            x, lo, hi = s["paper_ic"], s["paper_ic_lo"], s["paper_ic_hi"]
        else:
            s = d[d["background"] == bg].set_index("group").reindex(groups)
            x, lo, hi = s["ic"], s["ic025"], s["ic975"]
        ax.errorbar(x, y0 + o, xerr=[x - lo, hi - x], fmt=mk, color=col, ms=4.5, lw=1.2, capsize=2, label=label)
    ax.axvline(0, color="grey", lw=0.8, ls=":")
    ax.axvline(2, color="grey", lw=0.8, ls="--")
    ax.set_yticks(y0)
    ax.set_yticklabels(groups)
    ax.set_xlabel("Information Component (95% interval)")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.set_title("Fluoroquinolones and peripheral neuropathy: validation pull, rule B", fontsize=10)
    fig.tight_layout()
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "validation_forest.png", dpi=300)
    fig.savefig(out / "validation_forest.pdf")
    plt.close(fig)


def main():
    f = ROOT / "data" / "interim" / "validation" / "signals.csv"
    if not f.exists():
        sys.exit(f"Missing {f}. Run script 05 on the validation pull first.")
    t = pd.read_csv(f)
    tab, dec = build(t)
    tdir = ROOT / "results" / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    tab.to_csv(tdir / "validation_table.csv", index=False, float_format="%.4f")
    dec.to_csv(tdir / "validation_decomposition.csv", index=False, float_format="%.4f")
    forest(tab, ROOT / "results" / "figures")

    pd.set_option("display.width", 200)
    for rule in RULES:
        print(f"\nRule {rule}, background 'all': ours against Table 1")
        s = tab[(tab["rule"] == rule) & (tab["background"] == "all")]
        print(s[["group", "paper_O", "O", "paper_n_drug", "n_drug", "paper_ic", "ic", "d_ic"]].round(2).to_string(index=False))
    for rule in RULES:
        for bg in ("all", "no_concomitant"):
            print(f"\nRule {rule}, background '{bg}': gap decomposition (IC units)")
            s = dec[(dec["rule"] == rule) & (dec["background"] == bg)]
            print(s[["group", "paper_ic", "case_count_effect", "comparator_effect", "rounding", "total_gap"]]
                  .round(2).to_string(index=False))
    print(f"\nTables in {tdir}, figure in {ROOT / 'results' / 'figures'}")


if __name__ == "__main__":
    main()
