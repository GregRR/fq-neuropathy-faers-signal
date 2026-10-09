# Methods

Status: in progress. This file records what has been built and verified. Sections
for the signal-detection steps are added as each script is written.

## Aim

Replicate the BCPNN disproportionality analysis of Shamim et al. (2026) for
fluoroquinolone-associated peripheral neuropathy using public FAERS data, then run
the same pipeline on quarters the paper did not cover. This is a replication and
extension of published work. FAERS is spontaneous-report data, so the results
describe reporting patterns, not incidence or causation.

## Data source

FAERS quarterly ASCII extracts from the FDA, `$`-delimited text. Two analysis windows,
each processed separately:

| Pull | Quarters | Purpose |
|---|---|---|
| validation | 2007 Q1 to 2024 Q3 | Matches the window of Shamim et al. |
| extension | 2024 Q4 to 2026 Q2 | Data the paper did not see |

Download: `scripts/00_download_faers_quarters.py`.

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

## Case deduplication

Applied within each pull:

1. Remove any case listed in a Deleted file of a quarter in the pull.
2. Keep one report per CASEID. FAERS-era reports take priority over legacy reports.
   Within FAERS, keep the highest CASEVERSION (ties: latest FDA_DT, then highest
   PRIMARYID). Within legacy, keep the latest FDA_DT (ties: highest ISR).
3. Take DRUG and REAC rows from the same quarterly file as the retained report.

Each pull is deduplicated on its own, so the validation pull only contains data FDA
had released through 2024 Q3.

Differences from the published study: Shamim et al. used the FAERS Public Dashboard
and deduplicated on matching age, sex, country, initial FDA receipt date, drug,
indication and reaction (WHO Uppsala Monitoring Centre guidelines). The quarterly
extracts have no equivalent preprocessing, and legacy reports have no initial FDA
receipt date. Case counts will therefore differ from the paper's 37,778.

## Pending

Drug filtering, reaction filtering, contingency tables, BCPNN and ROR, validation
against the paper's Table 1, extension pull, sensitivity analyses.
