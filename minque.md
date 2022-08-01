---
title: MINQUE/REML
layout: template
filename: minque.md
---

## Input files
MPH uses the same file formats as [SSGP](https://sites.google.com/view/ssgp), another program that the same author developed. Refer to [this page](https://sites.google.com/view/ssgp/documentation/manual#h.p_QS3vj5saXQJZ) for details. In short, MPH prefers PLINK bim/fam/bed files for genotypes, and other input files are mostly CSV.

---

## Partitioning SNP heritability
```
mph --minque --grm_list chr.grms.txt --phenotype phen.csv --trait milk --num_threads 10 --output milk.chr
```

### To include covariates, add \-\-*covariate_file* and \-\-*covariate_names*.
```
--covariate_file covar.csv --covariate_names all
```
If *all* is specified, MPH will use as covariates all columns from the 2nd to the last in the covariate file.
```
--covariate_file covar.csv --covariate_names g1,g2,g3
```

### For DYD-like data, add \-\-*error_weight_name* to specify individual reliabilies.
```
--error_weight_name milk_wt
```
The error weights can be set to 1/*r*<sup>2</sup>-1.

### Other optional arguments
```--heritability 0.5```
The SNP heritability value for initializing MINQUE iterations. An accurate value may improve convergence. The default is 0.5.

```--num_iterations 20```
Max number of MINQUE iterations. The default is 20.

```--tol 0.01```
Tolerance. MINQUE iterations stop when logLL has a change smaller than that. The default is 0.01.

```--num_random_vectors 100```
Number of random vectors. The default (100) is usually sufficient.
