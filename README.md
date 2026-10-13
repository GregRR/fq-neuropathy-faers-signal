# fq-neuropathy-faers-signal

Independent replication and extension of a published FAERS pharmacovigilance
disproportionality analysis linking fluoroquinolone antibiotics (ciprofloxacin,
levofloxacin, moxifloxacin, ofloxacin, gemifloxacin, delafloxacin) to peripheral
neuropathy.

**Status:** the pipeline runs from download to the extension summary (scripts 00 to 07)
on both data pulls, and the control-drug and comparator-antibiotic analyses (scripts 08
and 09) have run on both. A comparison against all other systemic antibiotics is deferred. Validation against the published paper is partial: the estimator
reproduces the paper's Table 1 from the published counts, but counts and Information
Components rebuilt from raw FAERS files are lower than the published ones (see
Current results).

## What this is, and isn't

This project replicates a method, it does not claim a novel finding. A peer-reviewed
study (Shamim et al., *Clinical and Translational Science*, 2026) already applied a
Bayesian disproportionality analysis to FAERS data through 2024 Q3 and found a
significant signal between several fluoroquinolones and peripheral neuropathy-related
adverse event reports. The goal here is to demonstrate that method end to end, on real
public data, with documented and reproducible code, and to extend the published work in
three ways:

1. **Time extension.** Run the identical pipeline on FAERS data released after the
   paper's 2024 Q3 cutoff (2024 Q4 to 2026 Q2, new cases only), to see whether the
   signal holds, strengthens, or shifts.
2. **Method cross-check.** Compute a Reporting Odds Ratio (ROR) and a Bayesian
   Confidence Propagation Neural Network Information Component (IC), and compare this
   pipeline's per-drug results with the paper's Table 1, splitting any gap into the part
   due to case counts and the part due to the comparison population.
3. **Control drugs and comparator antibiotics.** Run azithromycin (negative control) and
   metronidazole (positive control) through the same pipeline, and compare each
   fluoroquinolone group with four comparator antibiotics (amoxicillin with or without
   clavulanate, cephalexin, sulfamethoxazole-trimethoprim, azithromycin) as well as with
   all of FAERS.

FAERS is spontaneous-reporting data. Per the FDA's own documentation, it cannot be used
to calculate incidence or to directly compare drug risk. Nothing in this repository
should be read as a claim of causation.

## Data source

FAERS ASCII quarterly extract files, public, free, no application or institutional
affiliation required. Available from the FDA at
https://www.fda.gov/drugs/fda-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-ascii-files.

FAERS data is U.S. government public domain. No separate data license applies. Raw and
intermediate data are not stored in this repository.

## Method

Details are in `docs/methods.md`. In brief:

- Cases: a fluoroquinolone as primary or secondary suspect, under two readings of "sole
  suspected drug" (rule A: no other suspect drug; rule B: no other drug in any role), with
  the paper's duplicate-removal rule applied in two variants.
- Event: 24 MedDRA preferred terms for peripheral neuropathy named in the paper's Figure 2
  (PN24), and the seven terms the paper retained (PN7).
- Comparison populations: all reports, reports without any fluoroquinolone, reports
  without any concomitant drug, and all reports with the event fraction fixed at the value
  implied by the paper's Table 1 (1.826%). The paper says only "the full FAERS dataset".
- Statistics: ROR with a Woolf interval and a Fisher exact p-value; IC from the Dirichlet
  posterior of Norén et al. (2006) with 2.5 and 97.5 percentile Monte Carlo limits. Two
  signal rules are reported side by side: IC025 above 0, and the paper's stated rule of
  IC above 2 with IC025 above 0.
- Control drugs: azithromycin and metronidazole, with the same case rules, duplicate
  removal, event and backgrounds.

## Current results

Validation pull (2007 Q1 to 2024 Q3), rule B, duplicate variant 1, all cases, PN24 event.
IC is the point estimate. "All reports" and "no concomitant" are the two comparison
populations defined above.

| Group | Paper (cases/reports, IC) | Here (cases/reports) | IC, all reports | IC, no concomitant |
|---|---|---|---|---|
| All fluoroquinolones | 5,588/37,778, 3.02 | 5,480/41,353 | 2.51 | 2.84 |
| Ciprofloxacin | 2,291/14,434, 3.12 | 2,310/15,732 | 2.66 | 2.99 |
| Levofloxacin | 1,942/13,598, 2.97 | 1,918/14,514 | 2.51 | 2.83 |
| Moxifloxacin | 642/7,439, 2.24 | 631/8,836 | 1.62 | 1.94 |
| Combination | 652/1,505, 4.54 | 569/1,452 | 4.06 | 4.38 |
| Ofloxacin | 46/588, 2.05 | 38/592 | 1.43 | 1.75 |
| Gemifloxacin | 11/142, 1.89 | 8/121 | 1.36 | 1.63 |
| Delafloxacin | 4/72, 1.31 | 6/106 | 1.13 | 1.40 |

- The IC calculation reproduces every IC in the paper's Table 1 to within 0.01 when given
  the published counts and a background event fraction of 1.826%.
- The PN24 share of all reports here is 2.32%, against 1.83% implied by the paper. That
  difference alone lowers every IC by 0.35. Reports without concomitant drugs have a
  PN24 share of 1.85%, close to the implied value.
- The remaining gap comes from case counts: here a smaller share of fluoroquinolone
  reports carries a neuropathy term (13.3% against 14.8% for all fluoroquinolones). The
  difference is largest for moxifloxacin (8,836 reports here, 7,439 in the paper).
- Group ordering matches the paper. Every group except delafloxacin has IC025 above 0 on
  the all-reports background. Moxifloxacin and ofloxacin fall below the paper's IC above 2
  cut-off. No source for that cut-off was found in the papers it cites.
- Extension pull (2024 Q4 to 2026 Q2): 304 neuropathy reports among 2,463 all-fluoroquinolone
  reports (12.3%), IC 2.62 against 2.51 in the validation pull, on the all-reports
  background. The background neuropathy share fell from 2.32% to 1.99%, so IC comparisons
  between pulls should be read alongside the neuropathy share among fluoroquinolone
  reports. Gemifloxacin has no extension reports, and delafloxacin has two.

- Control drugs, validation pull, rule B, on the all-reports background. Azithromycin
  (negative control): 192 neuropathy reports among 8,885, IC -0.10 (IC025 -0.31), ROR 0.93
  (0.80 to 1.07). Metronidazole (positive control): 545 among 5,732, IC 2.03 (IC025 1.91),
  ROR 4.42 (4.05 to 4.83). On the no-concomitant background azithromycin reaches IC 0.22
  (IC025 0.01, ROR 1.17), at the edge of the IC025 above 0 threshold and far below IC 2;
  metronidazole reaches 2.36. Metronidazole's IC falls between levofloxacin (2.51) and
  moxifloxacin (1.62). Extension pull: azithromycin 10/596, IC -0.23; metronidazole
  78/605, IC 2.65.
- Comparator antibiotics, validation pull, rule B. Neuropathy share among comparator
  reports: amoxicillin 1.3% (202/14,977), cephalexin 2.0% (40/1,972),
  sulfamethoxazole-trimethoprim 2.7% (170/6,397), azithromycin 2.2% (192/8,885), the four
  pooled 1.9% (604/32,231), against 13.3% for all fluoroquinolones. The ROR of all
  fluoroquinolones against the pooled comparators is 8.00 (7.34 to 8.71); by single
  comparator it runs from 5.60 (sulfamethoxazole-trimethoprim) to 11.17 (amoxicillin), so
  the size of the contrast depends on the comparator by about a factor of two. Every
  interval lies above 1 except delafloxacin against sulfamethoxazole-trimethoprim (2.20,
  0.95 to 5.08). Extension pull: ROR 9.85 (7.19 to 13.49) against the pooled comparators;
  levofloxacin is 7.69 against 7.97 in the validation pull, while the pooled comparator
  neuropathy share fell from 1.9% to 1.4%. Relative to the comparators the extension
  levofloxacin result does not drop, unlike its IC against all of FAERS.

Open questions about the gap with the paper are listed in `docs/methods.md`.

## Repository structure

```
scripts/
  00_download_faers_quarters.py          download the quarterly ZIPs
  00b_audit_faers_zips.py                check the downloads
  01_parse_and_link_core_files.py        parse DEMO, DRUG, REAC; remove repeated case versions
  02_filter_fq_and_comparator_drugs.py   flag the six fluoroquinolones; case rules A and B
  02b_deduplicate_fq_cases.py            paper-style duplicate removal
  03_filter_neuropathy_reactions.py      flag neuropathy terms; count them in all reports
  03b_flag_topical_cases.py              flag topical-only cases (sensitivity)
  04_build_contingency_tables.py         2x2 counts under four comparison populations
  05_compute_ror_and_bcpnn.py            ROR and IC with intervals; self-test against Table 1
  06_validate_against_table1.py          compare with the paper; decompose the gap
  07_summarise_extension.py              extension pull; pooled estimates; plots
  08_comparator_controls.py              control drugs and comparator antibiotics: case counts
  09_fq_vs_comparators.py                fluoroquinolone groups against comparator antibiotics
docs/
  methods.md
  references.md
results/
  tables/                                aggregate counts and estimates
  figures/
```

Each script's docstring gives its usage, inputs and outputs. Run from the repository
root, validation pull first, then the extension pull, for example:

```
python3 scripts/05_compute_ror_and_bcpnn.py validation
python3 scripts/05_compute_ror_and_bcpnn.py extension
python3 scripts/06_validate_against_table1.py
python3 scripts/07_summarise_extension.py
python3 scripts/08_comparator_controls.py validation
python3 scripts/08_comparator_controls.py extension
python3 scripts/08_comparator_controls.py summary
python3 scripts/09_fq_vs_comparators.py
```

The environment is in `environment.yml`.

## Benchmark papers

- Shamim N, Doughty K, Chau HT, Brown J, Baldock R, Au NPB. "Disproportionality
  Analysis of Fluoroquinolone-Associated Peripheral Neuropathy in the FAERS Database."
  *Clin Transl Sci.* 2026;19:e70541.
- Morales D, Pacurariu A, Slattery J, Pinheiro L, McGettigan P, Kurz X. "Association
  Between Peripheral Neuropathy and Exposure to Oral Fluoroquinolone or
  Amoxicillin-Clavulanate Therapy." *JAMA Neurol.* 2019;76(7):827-833.

Full list in `docs/references.md`.

## Companion project

This repository is the signal-detection half of a two-part body of work. The
companion repository, [fq-neuropathy-mimic-cohort](https://github.com/GregRR/fq-neuropathy-mimic-cohort),
answers this project's main stated limitation, FAERS data cannot support confounder
adjustment, with a propensity-score-matched cohort analysis in MIMIC-IV.

## License

MIT. See `LICENSE`.
