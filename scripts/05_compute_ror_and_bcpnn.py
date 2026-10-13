"""
05_compute_ror_and_bcpnn.py

Reporting odds ratio (ROR) and Bayesian Confidence Propagation Neural Network Information Component (IC)
for every row of contingency.csv (script 04).

Usage (from the repository root, after script 04):
    python3 scripts/05_compute_ror_and_bcpnn.py validation
    python3 scripts/05_compute_ror_and_bcpnn.py extension
    python3 scripts/05_compute_ror_and_bcpnn.py selftest

Information Component (Noren et al. 2006, Stat Med 25:3740-3757, doi:10.1002/sim.2473, Sections 2 to 2.3):
  counts n11 = O (drug and event), n1. = n_drug, n.1 = n_event, n.. = N
  prior Di(a11, a10, a01, a00) with a_ij = q_i. q_.j a.. and a.. = 0.5 / (q1. q.1)  (equations 3 to 5), so a11 = 0.5
  posterior Di(g11, g10, g01, g00), g_ij = a_ij + n_ij
  central estimate (equation 6): IC = log2( E[p11] / (E[p1.] E[p.1]) ), E[.] the posterior means
  interval: 2.5 and 97.5 percentiles of log2( p11 / (p1. p.1) ) over Monte Carlo draws from the posterior Dirichlet,
            the accuracy benchmark of the same paper (Sections 2 and 2.3)
Signal rules reported side by side:
  ic025_gt0   IC025 > 0, the routine screening threshold of the WHO database (Noren 2006, Section 2.5)
  paper_rule  IC > 2 and IC025 > 0, the rule stated by Shamim et al. 2026 (its cited sources, Chen et al. 2024 and
              Chau and Au 2026, state IC025 > 0 and an ROR rule respectively)
Reporting odds ratio: a*d / (b*c), 95% interval from the standard error of the log odds ratio (Woolf), with 0.5 added
  to every cell when any cell is zero. p value: two-sided Fisher exact test on the same 2x2 table.
  scipy returns exactly 0 below about 1e-308, so those p values are recomputed on the log scale from the
  hypergeometric distribution and stored as log10_p_fisher (always filled).

The selftest recomputes the rows of Table 1 of Shamim et al. 2026 using the event fraction
implied by the paper (1.826%) and prints them next to the published values.

Rows with no cases, or with a negative cell (an impossible table), get missing values.

Output in data/interim/<pull>/: signals.csv, the rows of contingency.csv plus ror, ror_lo, ror_hi, p_fisher, log10_p_fisher,
ic, ic025, ic975 and the signal flags.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import fisher_exact, hypergeom

ROOT = Path(__file__).resolve().parent.parent
N_DRAWS = 100_000
SEED = 20261009


def ic_posterior(n11, n1_, n_1, n, rng=None, draws=N_DRAWS):
    """Central IC (Noren 2006 eq 6) and the 2.5 / 97.5 percentiles from Monte Carlo draws."""
    n10, n01 = n1_ - n11, n_1 - n11
    n00 = n - n1_ - n_1 + n11
    q1, q0 = (n1_ + 0.5) / (n + 1), (n - n1_ + 0.5) / (n + 1)
    c1, c0 = (n_1 + 0.5) / (n + 1), (n - n_1 + 0.5) / (n + 1)
    a = 0.5 / (q1 * c1)
    g = np.array([q1 * c1 * a + n11, q1 * c0 * a + n10, q0 * c1 * a + n01, q0 * c0 * a + n00], dtype="float64")
    tot = g.sum()
    ic = np.log2((g[0] / tot) / (((g[0] + g[1]) / tot) * ((g[0] + g[2]) / tot)))
    if rng is None:
        return ic, np.nan, np.nan
    s = rng.gamma(g[:, None], size=(4, draws))
    p = s / s.sum(axis=0)
    d = np.log2(p[0] / ((p[0] + p[1]) * (p[0] + p[2])))
    lo, hi = np.percentile(d, [2.5, 97.5])
    return ic, lo, hi


def fisher_log10p(a, b, c, d):
    """Two-sided Fisher exact p value on the log10 scale (same rule as scipy: tables with probability
    no larger than the observed table, within a relative tolerance of 1e-7)."""
    a, b, c, d = int(a), int(b), int(c), int(d)
    m, n_event, n_drug = a + b + c + d, a + c, a + b
    k = np.arange(max(0, n_drug + n_event - m), min(n_drug, n_event) + 1)
    lp = hypergeom.logpmf(k, m, n_event, n_drug)
    obs = hypergeom.logpmf(a, m, n_event, n_drug)
    return min(0.0, logsumexp(lp[lp <= obs + np.log1p(1e-7)]) / np.log(10))


def fisher_p(a, b, c, d):
    """(p, log10 p). p is scipy's value; log10 p comes from the log-scale computation when p underflows."""
    p = fisher_exact([[a, b], [c, d]])[1]
    lg = fisher_log10p(a, b, c, d) if p < 1e-300 else float(np.log10(p))
    return p, lg


def fmt_p(lg):
    """Format a p value from its log10, e.g. 4.60e-16 or 3.10e-412."""
    e = int(np.floor(lg))
    return f"{10 ** (lg - e):.2f}e{e}"


def ror(a, b, c, d):
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    est = a * d / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return est, est * np.exp(-1.96 * se), est * np.exp(1.96 * se)


def main(pull):
    d = ROOT / "data" / "interim" / pull
    if not (d / "contingency.csv").exists():
        sys.exit(f"Missing {d / 'contingency.csv'}. Run script 04 first.")
    t = pd.read_csv(d / "contingency.csv")
    rng = np.random.default_rng(SEED)
    cols = {k: [] for k in ("ror", "ror_lo", "ror_hi", "p_fisher", "log10_p_fisher", "ic", "ic025", "ic975")}
    for i, r in enumerate(t.itertuples(index=False)):
        if r.O == 0 or r.n_drug == 0 or min(r.a, r.b, r.c, r.d) < 0:
            for k in cols:
                cols[k].append(np.nan)
            continue
        e = ror(r.a, r.b, r.c, r.d)
        ic = ic_posterior(r.O, r.n_drug, r.n_event, r.N, rng)
        p = fisher_p(r.a, r.b, r.c, r.d)
        for k, v in zip(cols, (*e, *p, *ic)):
            cols[k].append(v)
        if (i + 1) % 500 == 0:
            print(f"{i + 1:,} / {len(t):,} rows", flush=True)
    t = t.assign(**cols)
    t["ic025_gt0"] = t["ic025"] > 0
    t["paper_rule"] = (t["ic"] > 2) & (t["ic025"] > 0)
    t.to_csv(d / "signals.csv", index=False)

    sel = t[(t["rule"] == "B_no_other_drug_any_role") & (t["dedup"] == "missing_as_value") &
            (t["topical"] == "all_cases") & (t["event"] == "pn24")]
    print("\nRule B, duplicate variant 1, all cases, PN24 event:")
    for bg, g in sel.groupby("background"):
        print(f"\n  background {bg}")
        for r in g.itertuples():
            if pd.isna(r.ic):
                print(f"    {r.group:14s} {r.O:>6,}/{r.n_drug:<7,}  no estimate")
            else:
                print(f"    {r.group:14s} {r.O:>6,}/{r.n_drug:<7,} IC {r.ic:5.2f} ({r.ic025:5.2f}, {r.ic975:5.2f})  "
                      f"ROR {r.ror:5.2f} ({r.ror_lo:5.2f}, {r.ror_hi:5.2f})  p {fmt_p(r.log10_p_fisher)}")
    print(f"\n{len(t):,} rows written to {d / 'signals.csv'}")


def selftest():
    rng = np.random.default_rng(SEED)
    n = 17_610_899
    n_event = round(0.01826 * n)
    paper = {"all FQs": (5588, 37778, "3.02 (2.98, 3.05)", "<0.0001"),
             "ciprofloxacin": (2291, 14434, "3.12 (3.06, 3.17)", "<0.0001"),
             "levofloxacin": (1942, 13598, "2.97 (2.91, 3.02)", "<0.0001"),
             "moxifloxacin": (642, 7439, "2.24 (2.13, 2.34)", "<0.0001"),
             "combination": (652, 1505, "4.54 (4.46, 4.62)", "<0.0001"),
             "ofloxacin": (46, 588, "2.05 (1.62, 2.42)", "4.60E-16"),
             "gemifloxacin": (11, 142, "1.89 (0.95, 2.59)", "6.71E-05"),
             "delafloxacin": (4, 72, "1.31 (-0.39, 2.36)", "0.04")}
    print(f"Table 1 rows recomputed with N = {n:,} and an event fraction of 1.826% ({n_event:,} events)\n")
    for name, (o, nd, pub, pp) in paper.items():
        ic, lo, hi = ic_posterior(o, nd, n_event, n, rng, draws=1_000_000)
        p, lg = fisher_p(o, nd - o, n_event - o, n - nd - (n_event - o))
        print(f"  {name:13s} {o}/{nd:<6d} computed IC {ic:5.2f} ({lo:5.2f}, {hi:5.2f}) p {fmt_p(lg)}   published {pub}  p {pp}")


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) == 2 else ""
    if arg == "selftest":
        selftest()
    elif arg in ("validation", "extension"):
        main(arg)
    else:
        sys.exit("usage: 05_compute_ror_and_bcpnn.py [validation|extension|selftest]")
