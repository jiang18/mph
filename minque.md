---
title: MINQUE/REML
layout: template
filename: minque.md
---

## Input files
- \-\-grm_list: a space-delimited text file without header.
  - The first column lists GRM file path. 
  - The second column is optional and can list a label or an initial VC value for each GRM.
- \-\-phenotype: a CSV file with a header line. 
  - The first column must be the individual ID.
- \-\-covariate_file: a CSV file with a header line.
  - The first column must be the individual ID.
  - The covariate file is optional. 

Missing values of phenotypes or covariates should be left blank. Do not use -9, NA, or NaN. 

---

## REML or iterative MINQUE
```--minque```

```
mph --minque --grm_list chr.grms.txt --phenotype phen.csv --trait milk --num_threads 10 --output milk.chr
```

## MINQUE(0) or MINQUE(1)
To force MINQUE(0) or MINQUE(1), add *\-\-num_iterations 1*. The second column of the GRM list file should be set to 0 for MINQUE(0) and 1 for MINQUE(1).
```
mph --minque --num_iterations 1 --grm_list chr.grms.txt --phenotype phen.csv --trait milk --num_threads 10 --output milk.chr
```

## Memory-saving mode
```--save_memory```

---

## Covariates
To include covariates, add \-\-*covariate_file* and \-\-*covariate_names*.
```
--covariate_file covar.csv --covariate_names all
```
If *all* is specified, MPH will use as covariates all columns from the 2nd to the last in the covariate file.
```
--covariate_file covar.csv --covariate_names g1,g2,g3
```

## Deregressed estimated breeding values
For DYD-like data, add \-\-*error_weight_name* to specify individual reliabilies.
```
--error_weight_name milk_wt
```
The error weights can be set to 1/*r*<sup>2</sup>-1.

---

## General relationship matrix
The \-\-grm_list file can include any **general** relationship matrix, not just **genomic** relationship matrices.

---

## Additional options
```--heritability 0.5```
The SNP heritability value for initializing MINQUE iterations. An accurate value may improve convergence. The default is 0.5.

```--num_iterations 20```
Max number of MINQUE iterations. The default is 20.

```--tol 0.01```
Tolerance. MINQUE iterations stop when logLL has a change smaller than that. The default is 0.01.

```--num_random_vectors 100```
Number of random vectors. The default (100) is usually sufficient.
