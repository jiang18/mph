---
title: Download
layout: template
filename: download.md
---

## Software download
[https://github.com/jiang18/mph/releases/tag/20220715](https://github.com/jiang18/mph/releases/tag/20220715)

## Example data
[The QTL-MAS 2012 data](https://github.com/jiang18/mph/raw/main/QTL-MAS-2012.zip)  
Refer to [https://pubmed.ncbi.nlm.nih.gov/25519515/](https://pubmed.ncbi.nlm.nih.gov/25519515/) for how the data set was simulated.

---

## Update log
### July 15, 2022
- Added --deduct_grms for deducting GRMs.
- Added --dominance for computing dominance GRM.
- Added --make_core for computing covariance between random effects (CORE).
- Added --make_fore for computing first-order interaction between random effects (FORE).
- Improved --predict by outputting EBLUPs for each variance component.
- Improved --simulate by outputting the last 10% individuals' BLUPs.

### January 20, 2022
- Scaled logLL relative to null model.
- Added the sampling variance-covariance matrix of all VC estimates in the VC output file.

### September 6, 2021
- Fixed a bug in reading CSV files.

### September 2, 2021
- Added an output file to keep key information of --minque iterations.

### August 21, 2021
- Added an option '--'predict to predict total genetic values.
