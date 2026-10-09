# References

## Benchmark studies

- Shamim N, Doughty K, Chau HT, Brown J, Baldock R, Au NPB. Disproportionality
  Analysis of Fluoroquinolone-Associated Peripheral Neuropathy in the FAERS Database.
  *Clinical and Translational Science.* 2026;19:e70541. https://doi.org/10.1111/cts.70541
  Validation benchmark: BCPNN/IC per drug, FAERS Public Dashboard, 2007 to 2024 Q3.
- Morales D, Pacurariu A, Slattery J, Pinheiro L, McGettigan P, Kurz X. Association
  Between Peripheral Neuropathy and Exposure to Oral Fluoroquinolone or
  Amoxicillin-Clavulanate Therapy. *JAMA Neurology.* 2019;76(7):827-833.
  doi:10.1001/jamaneurol.2019.0887. Independent confirmation: UK THIN nested
  case-control, 1999 to 2015.

## Data and documentation

- FDA FAERS Quarterly Data Extract Files (download page):
  https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html
- FDA file layout documentation shipped inside each quarterly ZIP: `ASC_NTS.pdf`
  (data element definitions; section G maps legacy AERS fields to FAERS fields) and
  `Readme.pdf` (case versioning, the 2012 transition from legacy AERS to FAERS).

## Methods papers to retrieve before implementing BCPNN

Not yet read, so no details are recorded here.

- Bate et al. 1998 (original Bayesian Confidence Propagation Neural Network method)
- Noren et al. 2006 (shrinkage refinement of the Information Component)
