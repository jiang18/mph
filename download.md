---
title: Download
layout: template
filename: download.md
---

## Software download
[https://github.com/jiang18/mph/releases/tag/20220908](https://github.com/jiang18/mph/releases/tag/20220908)

## Example data
[The QTL-MAS 2012 data](https://github.com/jiang18/mph/raw/main/QTL-MAS-2012.zip)  
Refer to [https://pubmed.ncbi.nlm.nih.gov/25519515/](https://pubmed.ncbi.nlm.nih.gov/25519515/) for how the data set was simulated.

---

## Update log
### Sep 08, 2022
- Added a memory-saving mode (\-\-save_ram).
- Improved GRM reading by multithreading.

### Aug 22, 2022
- Enabled on-the-fly genotype reading for making GRM.

### Jul 29, 2022
- Improved the computation of MAFs and HWE p-values by multithreading.

### Jul 15, 2022
- Added \-\-deduct_grms for deducting GRMs.
- Added \-\-dominance for computing dominance GRM.
- Added \-\-make_core for computing COvariance between Random Effects (CORE).
- Added \-\-make_fore for computing First-Order interaction between Random Effects (FORE).
- Improved \-\-predict by outputting EBLUPs for each variance component.
- Improved \-\-simulate by outputting the last 10% individuals' BLUPs.

### Jan 20, 2022
- Scaled logLL relative to null model.
- Added the sampling variance-covariance matrix of all VC estimates in the VC output file.

### Sep 06, 2021
- Fixed a bug in reading CSV files.

### Sep 02, 2021
- Added an output file to keep key information of \-\-minque iterations.

### Aug 21, 2021
- Added an option \-\-predict to predict total genetic values.
