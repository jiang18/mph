---
title: MINQUE/REML
layout: template
filename: minque.md
---

## Table of contents
- [Input](#input)
- [Options](#options)
    - [REML or iterative MINQUE](#reml-or-iterative-minque)
    - [MINQUE(0) or MINQUE(1)](#minque0-or-minque1)
    - [Memory-saving mode](#memory-saving-mode)
    - [Covariates](#covariates)
    - [Deregressed estimated breeding values](#deregressed-estimated-breeding-values)
    - [General relationship matrix](#general-relationship-matrix)
    - [Additional options](#additional-options)
- [Output](#output)

## Input
* \-\-grm_list: a space-delimited text file without header.
  * The first column lists GRM file path. 
  * The second column is optional and can list an initial VC value for each GRM.
* \-\-phenotype: a CSV file with a header line. 
  * The first column must be the individual ID.
* \-\-covariate_file: a CSV file with a header line.
  * The first column must be the individual ID.
  * The covariate file is optional. 

Missing values of phenotypes or covariates should be left blank. Do not use -9, NA, or NaN. 

---

## Options

### REML or iterative MINQUE
```
--minque
```
This flag turns on iterative MINQUE (equivalent to the Fisher scoring algorithm for REML).

```
mph --minque --grm_list chr.grms.txt --phenotype phen.csv --trait milk --num_threads 10 --output milk.chr
```

### MINQUE(0) or MINQUE(1)
To force MINQUE(0) or MINQUE(1), add *\-\-num_iterations 1*. The second column of the GRM list file should be set to 0 for MINQUE(0) and 1 for MINQUE(1).
```
mph --minque --num_iterations 1 --grm_list chr.grms.txt --phenotype phen.csv --trait milk --num_threads 10 --output milk.chr
```

### Memory-saving mode
```
--save_memory
```
This flag enables the memory-saving mode. This mode is not necessarily slower, particularly on a solid-state drive.

---

### Covariates
To include covariates, add \-\-*covariate_file* and \-\-*covariate_names*.
```
--covariate_file covar.csv --covariate_names all
```
If *all* is specified, MPH will use as covariates all columns from the 2nd to the last in the covariate file.
```
--covariate_file covar.csv --covariate_names g1,g2,g3
```

### Deregressed estimated breeding values
For DYD-like data, add \-\-*error_weight_name* to specify individual reliabilies.
```
--error_weight_name milk_wt
```
The error weights can be set to 1/*r*<sup>2</sup>-1.

---

### General relationship matrix
The \-\-grm_list file can include any **general** relationship matrix, not just **genomic** relationship matrices.

---

### Additional options
```--heritability 0.5```
The SNP heritability value for initializing MINQUE iterations. An accurate value may improve convergence. The default is 0.5.

```--num_iterations 20```
Max number of MINQUE iterations. The default is 20.

```--tol 0.01```
Tolerance. MINQUE iterations stop when logLL has a change smaller than that. The default is 0.01.

```--num_random_vectors 100```
Number of random probing vectors for stochastic trace estimation. The default (100) is usually sufficient.

```--seed 0```
Random seed in stochastic trace estimation. The default is 0.

---

## Output
* mq.blue.csv: best linear unbiased estimates of fixed effects
* mq.iter.csv: summary of REML iterations
* mq.vc.csv: variance component estimates
| Column | Description |
|----------|----------|
| name | Name of the variance component or GRM |
| m | Weighted number of SNPs in the GRM |
| var | Estimate of the variance component |
| seV | Standard error of the variance component estimate |
| pve | Estimate of the proportion of variance explained |
| seP | Standard error of the estimate of the proportion of variance explained |
| enrichment | Estimate of the per-SNP heritability enrichment |
| seE | Standard error of the estimate of the per-SNP heritability enrichment |

* mq.py.csv: residuals of the linear mixed model
