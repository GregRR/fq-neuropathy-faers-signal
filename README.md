# fq-neuropathy-faers-signal

Independent replication and extension of a published FAERS pharmacovigilance
disproportionality analysis linking fluoroquinolone antibiotics (ciprofloxacin,
levofloxacin, moxifloxacin, ofloxacin, gemifloxacin, delafloxacin) to peripheral
neuropathy.

**Status: planning complete, pipeline build not yet started.**

## What this is, and isn't

This project replicates a method, it does not claim a novel finding. A peer-reviewed
study (Shamim et al., *Clinical and Translational Science*, 2026) already applied a
Bayesian disproportionality analysis to FAERS data through 2024 Q3 and found a
significant signal between several fluoroquinolones and peripheral neuropathy-related
adverse event reports. The goal here is to demonstrate that method end to end, on real
public data, with fully documented and reproducible code, and to extend the published
work in two specific, honest ways:

1. **Time extension.** Re-run the identical pipeline on FAERS data released after the
   published paper's 2024 Q3 cutoff, to see whether the signal holds, strengthens, or
   shifts.
2. **Method cross-check.** Compute both a classical Reporting Odds Ratio (ROR) and a
   from-scratch Bayesian Confidence Propagation Neural Network (BCPNN) implementation
   (the method the benchmark paper actually used), then validate this pipeline's
   per-drug results directly against the paper's published numbers.

FAERS is spontaneous-reporting data. Per the FDA's own documentation, it cannot be used
to calculate incidence or to directly compare drug risk. Nothing in this repository
should be read as a claim of causation.

## Data source

FAERS ASCII quarterly extract files, public, free, no application or institutional
affiliation required. Available from the FDA at
https://www.fda.gov/drugs/fda-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-ascii-files.

FAERS data is U.S. government public domain. No separate data license applies.

## Benchmark papers

- Shamim N, Doughty K, Chau HT, Brown J, Baldock R, Au NPB. "Disproportionality
  Analysis of Fluoroquinolone-Associated Peripheral Neuropathy in the FAERS Database."
  *Clin Transl Sci.* 2026;19:e70541.
- Morales D, Pacurariu A, Slattery J, Pinheiro L, McGettigan P, Kurz X. "Association
  Between Peripheral Neuropathy and Exposure to Oral Fluoroquinolone or
  Amoxicillin-Clavulanate Therapy." *JAMA Neurol.* 2019;76(7):827-833.

## Method

- Reporting Odds Ratio (ROR), computed from drug-event 2x2 contingency tables.
- Bayesian Confidence Propagation Neural Network (BCPNN), producing an Information
  Component (IC) per drug-event pair, following Bate et al. (1998) and Norén et al.
  (2006).
- Adverse event definition: seven MedDRA preferred terms confirmed from the benchmark
  paper (paraesthesia, peripheral neuropathy, hypoaesthesia, burning sensation,
  sensory disturbance, small fiber neuropathy, hyperaesthesia).
- Comparator drug set: a primary four-drug comparator (amoxicillin/amoxicillin-
  clavulanate, cephalexin, azithromycin, trimethoprim-sulfamethoxazole), with an
  all-other-antibiotics comparison as a secondary sensitivity check.

## Repository structure (planned)

```
00_download_faers_quarters.py
01_parse_and_link_core_files.py
02_filter_fq_and_comparator_drugs.py
03_filter_neuropathy_reactions.py
04_build_contingency_tables.py
05_compute_ror_and_bcpnn.py
06_validate_against_published_paper.py
07_extend_to_recent_quarters.py
08_sensitivity_vs_all_antibiotics.py
METHODS_NOTES.md
```

## Companion project

This repository is the signal-detection half of a two-part body of work. The
companion repository, [fq-neuropathy-mimic-cohort](https://github.com/GregRR/fq-neuropathy-mimic-cohort),
answers this project's main stated limitation, FAERS data cannot support confounder
adjustment, with a propensity-score-matched cohort analysis in MIMIC-IV.

## License

MIT. See `LICENSE`.
