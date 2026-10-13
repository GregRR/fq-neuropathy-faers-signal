"""
07_summarise_extension.py

Summarises the extension pull (new cases reported after the validation window, script 05) and sets it beside the
validation pull.

Usage (from the repository root, after script 05 on both pulls):
    python3 scripts/07_summarise_extension.py

Selection from signals.csv: duplicate variant 1 (missing_as_value), all cases, PN24 event, both sole-suspect rules.
The extension holds only cases not in the validation pull, so the two pulls are disjoint and their counts add.
The pooled rows add the counts of the two pulls (reports with the drug and the event, reports with the drug,
reports with the event in the background, reports in the background) and recompute ROR and IC on the sums.

Backgrounds: all, excl_fq_any and no_concomitant. The paper-implied background fixes the event fraction at the
value implied by Table 1 of Shamim et al. 2026 and has no meaning for a later pull, so it is not used here.

Output:
  results/tables/extension_summary.csv   one row per rule, background, group, pull (validation, extension, pooled)
  results/figures/extension_forest.png, .pdf
        IC with 95% interval, validation against extension, rule B, background 'all' (left) and
        'no_concomitant' (right). Groups with no estimate in the extension are marked n/a.
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
BACKGROUNDS = ("all", "excl_fq_any", "no_concomitant")
GROUPS = ("all", "ciprofloxacin", "levofloxacin", "moxifloxacin", "combination", "ofloxacin", "gemifloxacin",
          "delafloxacin")
KEEP = ["O", "n_drug", "n_event", "N", "ror", "ror_lo", "ror_hi", "ic", "ic025", "ic975"]


def pick(t, rule, bg):
    s = t[(t["rule"] == RULES[rule]) & (t["dedup"] == "missing_as_value") & (t["topical"] == "all_cases") &
          (t["event"] == "pn24") & (t["background"] == bg)]
    return s.set_index("group").reindex(list(GROUPS))


def pooled(v, e, rng):
    """Add the counts of two disjoint pulls and recompute ROR and IC. Rows with no drug reports in a pull
    contribute zero counts (the group is then the other pull's)."""
    o, nd = v["O"].fillna(0) + e["O"].fillna(0), v["n_drug"].fillna(0) + e["n_drug"].fillna(0)
    ne, n = v["n_event"] + e["n_event"], v["N"] + e["N"]
    rows = []
    for g in GROUPS:
        if o[g] == 0 or nd[g] == 0:
            rows.append(dict.fromkeys(KEEP, np.nan))
            continue
        a, b, c, d = o[g], nd[g] - o[g], ne[g] - o[g], n[g] - nd[g] - (ne[g] - o[g])
        r = s05.ror(a, b, c, d)
        ic = s05.ic_posterior(int(o[g]), int(nd[g]), int(ne[g]), int(n[g]), rng)
        rows.append(dict(zip(KEEP, (o[g], nd[g], ne[g], n[g], *r, *ic))))
    return pd.DataFrame(rows, index=list(GROUPS))


def build(tv, te):
    rng = np.random.default_rng(s05.SEED)
    out = []
    for rule in RULES:
        for bg in BACKGROUNDS:
            v, e = pick(tv, rule, bg), pick(te, rule, bg)
            p = pooled(v, e, rng)
            for pull, s in (("validation", v), ("extension", e), ("pooled", p)):
                s = s.reindex(columns=KEEP)
                for g in GROUPS:
                    out.append({"rule": rule, "background": bg, "group": g, "pull": pull, **s.loc[g].to_dict()})
    res = pd.DataFrame(out)
    res["event_fraction"] = res["n_event"] / res["N"]
    return res


def forest(res, out):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5), sharey=True)
    y0 = np.arange(len(GROUPS))[::-1].astype(float)
    series = (("Validation pull", "validation", "tab:blue", "s", 0.12),
              ("Extension pull", "extension", "tab:orange", "^", -0.12))
    for ax, bg, title in zip(axes, ("all", "no_concomitant"), ("Background 'all'", "Background 'no_concomitant'")):
        d = res[(res["rule"] == "B") & (res["background"] == bg)]
        for label, pull, col, mk, o in series:
            s = d[d["pull"] == pull].set_index("group").reindex(list(GROUPS))
            ok = s["ic"].notna().to_numpy()
            ax.errorbar(s["ic"][ok], y0[ok] + o, xerr=[(s["ic"] - s["ic025"])[ok], (s["ic975"] - s["ic"])[ok]],
                        fmt=mk, color=col, ms=4.5, lw=1.2, capsize=2, label=label)
            if pull == "extension":
                for yy in y0[~ok]:
                    ax.text(0.05, yy + o, "n/a", color=col, va="center", fontsize=8)
        ax.axvline(0, color="grey", lw=0.8, ls=":")
        ax.axvline(2, color="grey", lw=0.8, ls="--")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Information Component (95% interval)")
    axes[0].set_yticks(y0)
    axes[0].set_yticklabels(GROUPS)
    axes[0].legend(loc="lower right", fontsize=8, frameon=False)
    fig.suptitle("Fluoroquinolones and peripheral neuropathy: validation and extension pulls, rule B", fontsize=10)
    fig.tight_layout()
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "extension_forest.png", dpi=300)
    fig.savefig(out / "extension_forest.pdf")
    plt.close(fig)


def main():
    paths = {p: ROOT / "data" / "interim" / p / "signals.csv" for p in ("validation", "extension")}
    for f in paths.values():
        if not f.exists():
            sys.exit(f"Missing {f}. Run script 05 on both pulls first.")
    tv, te = (pd.read_csv(f) for f in paths.values())
    res = build(tv, te)
    tdir = ROOT / "results" / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    res.to_csv(tdir / "extension_summary.csv", index=False, float_format="%.4f")
    forest(res, ROOT / "results" / "figures")

    pd.set_option("display.width", 200)
    for rule in RULES:
        for bg in ("all", "no_concomitant"):
            print(f"\nRule {rule}, background '{bg}'")
            s = res[(res["rule"] == rule) & (res["background"] == bg)]
            w = s.pivot(index="group", columns="pull", values=["O", "n_drug", "ic", "ic025"]).reindex(list(GROUPS))
            w = w.reindex(columns=[(m, p) for m in ("O", "n_drug", "ic", "ic025")
                                   for p in ("validation", "extension", "pooled")])
            w.columns = [f"{m}_{p[:3]}" for m, p in w.columns]
            print(w.round(2).to_string())
    print("\nEvent fraction (PN24 reports / reports in the background), rule B:")
    f = res[(res["rule"] == "B") & (res["group"] == "all") & (res["pull"] != "pooled")]
    print(f.pivot(index="background", columns="pull", values="event_fraction").round(4).to_string())
    print(f"\nTable in {tdir}, figure in {ROOT / 'results' / 'figures'}")


if __name__ == "__main__":
    main()
