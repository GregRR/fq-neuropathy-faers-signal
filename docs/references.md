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
- Guo B, Zou Z, Lu J, Tan S, Mai L. Real-world safety signals of fluoroquinolone-associated
  peripheral neuropathy: analysis of US and Canadian databases combined with mechanistic
  exploration. *International Journal of Antimicrobial Agents.* 2026;67:107993.
  doi:10.1016/j.ijantimicag.2026.107993. Source of the azithromycin and metronidazole
  control drugs.

## Methods

- Bate A, Lindquist M, Edwards IR, Olsson S, Orre R, Lansner A, De Freitas RM. A Bayesian
  neural network method for adverse drug reaction signal generation. *European Journal of
  Clinical Pharmacology.* 1998;54:315-321. doi:10.1007/s002280050466.
- Norén GN, Bate A, Orre R, Edwards IR. Extending the methods used to screen the WHO drug
  safety database towards analysis of complex associations and improved accuracy for rare
  events. *Statistics in Medicine.* 2006;25(21):3740-3757. doi:10.1002/sim.2473.
  Prior, point estimate and interval used for the IC (equations 3 to 6, Sections 2 to 2.3
  and 2.5).
- Norén GN, Hopstadius J, Bate A. Shrinkage observed-to-expected ratios for robust and
  transparent large-scale pattern discovery. *Statistical Methods in Medical Research.*
  2013;22:57-69. Source of the shrinkage form IC = log2((O + 0.5) / (E + 0.5)).
- Chen Y, Ren X, Dai Y, Wang Y. Pharmacovigilance study of the association between
  peripheral neuropathy and antibody-drug conjugates using the FDA adverse event reporting
  system. *Scientific Reports.* 2024;14:21386. doi:10.1038/s41598-024-71977-0. Cited by
  Shamim et al. for the BCPNN method; states IC025 above 0 as the signal criterion.
- Chau HT, Au NPB. Fluoroquinolone-associated psychiatric and ocular adverse events: a
  disproportionality analysis using real-world data from FAERS (2011 to 2024).
  *Pharmacology Research & Perspectives.* 2026;14(1):e70206. doi:10.1002/prp2.70206.
  Cited by Shamim et al.; reporting odds ratio analysis only.

## Data and documentation

- FDA FAERS Quarterly Data Extract Files (download page):
  https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html
- FDA file layout documentation shipped inside each quarterly ZIP: `ASC_NTS.pdf`
  (data element definitions; section G maps legacy AERS fields to FAERS fields) and
  `Readme.pdf` (case versioning, the 2012 transition from legacy AERS to FAERS).
- Brand names for the control drugs and comparator antibiotics: drugs.com pages for
  azithromycin, metronidazole, topical metronidazole, vaginal metronidazole, Pylera,
  amoxicillin, amoxicillin and clavulanate, cephalexin and sulfamethoxazole and trimethoprim,
  and the drugs.com H. pylori eradication agents class page for multi-agent products
  (https://www.drugs.com).
  Brand names for the fluoroquinolones: drugs.com quinolone class page; topical products:
  Missouri HealthNet ophthalmic fluoroquinolone review and Nevada Medicaid otic
  fluoroquinolone review.
