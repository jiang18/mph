# Change Log

## July 16, 2024 (v0.54.0)
- Implemented custom genotype coding for GRM computation (--snp_genotype_coding).

## May 16, 2024 (v0.53.2)
- Added the sampling variance-covariance matrix for all fixed-effect estimates to the BLUE output file.
- Corrected the fixed-effect estimate of the intercept.

## Jan 29, 2024 (v0.52.1)
- Fixed a bug in --pred. 

## Jan 27, 2024 (v0.52.0)
- Updated the structure of the .mq.vc.csv file format to improve readability.
- Eliminated subroutines related to eigen-decomposition-based single-GRM REML.

## Nov 26, 2023 (v0.51.0)
- Implemented multi-trait REML.
- Fine-tuned the threshold for rank determination.

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
- Added the sampling variance-covariance matrix of all VC estimates to the VC output file.

## Sep 06, 2021
- Fixed a bug in reading CSV files.

## Sep 02, 2021
- Added an output file to keep key information of --minque iterations.

## Aug 21, 2021
- Added an option --predict to predict total genetic values.
