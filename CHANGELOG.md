# Change Log

## Sep 18, 2023 (v0.49.1)
- Added a flag (--exact_init) to initialize exact values for variance components.

## Apr 10, 2023 (v0.49.0)
- Enabled user-specified SNP weighting for dominance GRM.

## Sep 19, 2022
- Improved performance by reducing overhead.

## Sep 08, 2022
- Added a memory-saving mode (--save_memory).
- Improved GRM reading by multithreading.

## Aug 22, 2022
- Enabled on-the-fly genotype reading for making GRM.

## Jul 29, 2022
- Improved the computation of MAFs and HWE p-values by multithreading.

## Jul 15, 2022
- Added --deduct_grms for deducting GRMs.
- Added --dominance for computing dominance GRM.
- Added --make_core for computing COvariance between Random Effects (CORE).
- Added --make_fore for computing First-Order interaction between Random Effects (FORE).
- Improved --predict by outputting EBLUPs for each variance component.
- Improved --simulate by outputting the last 10% individuals' BLUPs.

## Jan 20, 2022
- Scaled logLL relative to null model.
- Added the sampling variance-covariance matrix of all VC estimates in the VC output file.

## Sep 06, 2021
- Fixed a bug in reading CSV files.

## Sep 02, 2021
- Added an output file to keep key information of --minque iterations.

## Aug 21, 2021
- Added an option --predict to predict total genetic values.