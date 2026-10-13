# Methods

Status: scripts 00 to 08 have been run on both pulls. The comparison with other
antibiotics as a background is not yet written.

## Aim

Replicate the BCPNN disproportionality analysis of Shamim et al. (2026) for
fluoroquinolone-associated peripheral neuropathy using public FAERS data, then run the
same pipeline on quarters the paper did not cover. This is a replication and extension of
published work. FAERS is spontaneous-report data, so the results describe reporting
patterns, not incidence or causation.

## Data source

FAERS quarterly ASCII extracts from the FDA, `$`-delimited text. Two analysis windows,
each processed separately:

| Pull | Quarters | Purpose |
|---|---|---|
| validation | 2007 Q1 to 2024 Q3 | Matches the window of Shamim et al. |
| extension | 2024 Q4 to 2026 Q2 | Data the paper did not see; cases first seen after 2024 Q3 only |

Download: `scripts/00_download_faers_quarters.py`; audit: `scripts/00b_audit_faers_zips.py`.

## File layouts

The extracts use two layouts, verified against the FDA's `ASC_NTS.pdf` and the header
rows of the files.

- Legacy AERS (through 2012 Q3): reports keyed on ISR; a case (CASE) can have
  several ISRs. No initial FDA receipt date and no active-ingredient field.
- FAERS (2012 Q4 onward): reports keyed on PRIMARYID (CASEID plus CASEVERSION).
  Sex column is named GNDR_COD until 2014 Q3 and SEX after. Active ingredient
  (PROD_AI) is available from 2014 Q3.

`scripts/01_parse_and_link_core_files.py` reads both layouts into one schema. Report
identifiers are prefixed `L` (legacy) or `F` (FAERS) so the two ID spaces cannot collide.
Because some legacy ISRs recur across quarters, cases are keyed on report identifier and
quarter together.

## Repeated case versions

Applied within each pull (script 01):

1. Remove any case listed in a Deleted file of a quarter in the pull.
2. Keep one report per CASEID. FAERS-era reports take priority over legacy reports.
   Within FAERS, keep the highest CASEVERSION (ties: latest FDA_DT, then highest
   PRIMARYID). Within legacy, keep the latest FDA_DT (ties: highest ISR).
3. Take DRUG and REAC rows from the same quarterly file as the retained report.

Each pull is deduplicated on its own, so the validation pull only contains data FDA had
released through 2024 Q3. The extension pull then drops every CASEID already present in
the validation pull.

## Fluoroquinolone identification (script 02)

DRUGNAME and PROD_AI are searched for generic-name stems and the U.S. brand names Cipro,
Proquin, Levaquin, Avelox, Floxin, Factive and Baxdela. Legacy files have no PROD_AI, so
a lookup is built from rows where the FDA's PROD_AI identifies an ingredient but DRUGNAME
matched none of the patterns (for example CRAVIT, TAVANIC, CIPROXIN). A name enters the
lookup only if it appears at least three times and PROD_AI maps it to one ingredient. One
name, MYLOTARG (gemtuzumab), carries a levofloxacin PROD_AI in five rows through a
source-data error and is excluded.

## Case definitions (scripts 02, 02b)

The paper says it kept reports where a fluoroquinolone was the sole suspected drug and
excluded reports with concomitant non-fluoroquinolone drugs, without saying which drug
roles it counted. Two readings are run side by side.

- Rule A: at least one fluoroquinolone as primary or secondary suspect, and no
  non-fluoroquinolone drug as primary or secondary suspect (concomitant drugs allowed).
- Rule B: at least one fluoroquinolone as primary or secondary suspect, and no
  non-fluoroquinolone drug in any role. Rule B is the primary reading.

Groups: all, each of the six ingredients (exactly one ingredient as suspect), and
combination (two or more).

Second-stage duplicate removal (script 02b) follows the paper: reports with the same age,
sex, country, initial FDA receipt date, drug, indication and reaction are duplicates. The
paper does not say how missing fields or the drug field were handled, so choices were made:

- Age is converted to years and rounded to 0.1; country is the occurrence country when
  present, else the reporter country; legacy reports use FDA_DT because they have no
  initial receipt date.
- Drug: the set of suspect fluoroquinolone ingredients, plus the set of every other drug
  on the report in any role (PROD_AI, else drug name). Indication: the set of indication
  terms across all drugs. Reaction: the set of reaction terms.
- Variant 1 (`missing_as_value`, primary): a missing age, sex or country counts as a
  value, so two reports missing the same field can match. Variant 2 (`conservative`):
  reports missing any of the three are never merged.
- Within a group of duplicates the lowest CASEID is kept. Duplicate removal is applied to
  fluoroquinolone cases only; the comparison populations are not deduplicated this way.

Topical products (script 03b): cases whose suspect fluoroquinolone rows are all topical
(eye, ear or skin products, identified by route or name) are flagged. The paper does not
say whether they were included. The scenario `excl_topical_only` drops them; it is a
sensitivity analysis.

## Event definition (script 03)

MedDRA preferred terms, normalised to UK spelling (hypoesthesia matches hypoaesthesia;
fiber matches fibre).

- PN24: the 24 terms named in Figure 2 of the paper (terms with at least three reports).
- PN7: the seven terms the paper retained as signals in Figure 2D (paraesthesia,
  neuropathy peripheral, hypoaesthesia, burning sensation, sensory disturbance, small
  fibre neuropathy, hyperaesthesia).

The paper selected 28 terms (its Table S9). Four had fewer than three reports and are not
named in the main text, so they are missing here. Table 1 of the paper counts cases with
any of the 28, so PN24 only approximates it.

## Comparison populations (script 04)

The paper says only "the full FAERS dataset". Four definitions are computed; all contain
the group being tested, adding it back where needed.

| Background | Definition |
|---|---|
| all | every retained report in the pull |
| excl_fq_any | every report with no fluoroquinolone in any role |
| no_concomitant | every report with no drug row in the concomitant role |
| paper_implied | all reports, with the event fraction fixed at 1.826% (PN24 only), the value that reproduces Table 1 |

Individual terms use the `all` background only.

## Statistics (script 05)

- ROR: a*d / (b*c) from the 2x2 table, 95% interval from the standard error of the log
  odds ratio (Woolf), 0.5 added to every cell when any cell is zero. Two-sided Fisher
  exact p-value, recomputed on the log scale when the value underflows (stored as
  `log10_p_fisher`).
- IC: Dirichlet posterior following Norén et al. (2006). Prior alpha_ij = q_i q_j alpha..
  with alpha.. = 0.5 / (q1. q.1) (equations 3 to 5), so alpha_11 = 0.5. The point estimate
  is log2(E[p11] / (E[p1.] E[p.1])) (equation 6). The 2.5 and 97.5 percentiles come from
  100,000 Monte Carlo draws from the posterior, the accuracy benchmark of the same paper,
  with a fixed seed.
- Self-test: with the paper's published counts and a background event fraction of 1.826%,
  the computed IC and its limits match Table 1 within 0.01 for all eight groups.
- Signal rules, both reported: IC025 above 0 (the routine screening threshold of the WHO
  database, Norén 2006, Section 2.5), and IC above 2 with IC025 above 0 (the rule stated
  by Shamim et al.). The shrinkage point estimate IC = log2((O + 0.5) / (E + 0.5)) of
  Norén, Hopstadius and Bate (2013), as used by Chen et al. (2024), which the paper cites
  for its method, gives similar values for large counts. No source for the requirement that IC exceed 2 was found in the papers the
  paper cites for it (Chen et al. 2024; Chau and Au 2026), nor in Norén et al. (2006).

## Validation against Table 1 (script 06)

The difference between the IC computed here and the published IC is split into:

- case-count effect: IC from our counts on the paper-implied background, minus IC from
  the published counts on the same background;
- comparator effect: IC from our counts on our background, minus IC from our counts on
  the paper-implied background;
- rounding: IC from the published counts recomputed here, minus the printed value.

The comparator effect is the same for every group, because it comes from the difference
in background event fraction (2.32% here on the `all` background against 1.83% implied).
On the `no_concomitant` background (1.85%) it nearly disappears.

## Extension (script 07)

The extension holds only cases not in the validation pull, so the pulls are disjoint and
their counts add. Pooled rows add the counts and recompute ROR and IC. The
`paper_implied` background has no meaning for a later pull and is not used. Comparisons
between pulls should look at the neuropathy share among fluoroquinolone reports as well
as the IC, since the background event fraction also changes between pulls.

## Control drugs (script 08)

Azithromycin (negative control) and metronidazole (positive control), following Guo et
al. (2026), who did not publish their drug definitions in the main text. Case rules,
duplicate removal, event and backgrounds are the same as for the fluoroquinolones, with
"other drug" meaning any drug outside the control group. Rule A is reported on the `all`
and `paper_implied` backgrounds only. Names are matched on generic stems, U.S. brand names
and a PROD_AI lookup as for the fluoroquinolones; the multi-ingredient product Pylera
(bismuth subcitrate potassium, metronidazole, tetracycline) is excluded from the
metronidazole group. Topical and vaginal products are flagged and handled as in script
03b.

## Differences from the published study

- The paper used the FAERS Public Dashboard and WHO Uppsala Monitoring Centre duplicate
  guidance. The quarterly extracts have no equivalent preprocessing, so duplicate removal
  here is a reconstruction from the paper's description.
- The paper's report counts (for example 37,778 for all fluoroquinolones under its
  case definition) are not reproduced exactly; rule B gives 41,353 in the validation pull.
- The event definition uses 24 of the paper's 28 terms.
- Comparison population, drug-role handling and the IC rule are not fully specified in the
  paper; see above.

## Open questions

- Which drug roles and which comparison population the paper used, and whether its
  denominator excluded reports with concomitant drugs. The implied event fraction
  (1.826%) is close to the `no_concomitant` fraction here (1.85%); the paper's total
  report count is not known.
- The source of the moxifloxacin excess (8,836 reports against 7,439; a candidate is the
  large number of U.S. single-agent Avelox reports from 2008, many with missing age), and
  the differences for gemifloxacin (121 against 142 reports), delafloxacin (106 against
  72), burning sensation and the combination group.
- Drift in MedDRA coding over the study window.
- The source of the IC above 2 requirement.

## Pending

Fluoroquinolone groups against a background of other systemic antibiotics.
