---
title: MINQUE
layout: template
filename: minque.md
---

## Input file formats
MPH uses the same file formats as [SSGP](https://sites.google.com/view/ssgp), another program that the same author developed. Refer to [this page](https://sites.google.com/view/ssgp/documentation/manual#h.p_QS3vj5saXQJZ) for details. In short, MPH prefers PLINK bim/fam/bed files for genotypes, and other input files are mostly CSV.

## SNP info file
This CSV file is critical for building genomic relationship matrices (GRMs). Each column corresponds to one GRM and each row corresponds to a SNP. In the example data, there are five columns for five chromosomes.

If a SNP does not belong to a GRM, leave the corresponding cell blank. If a SNP belongs to ***n*** GRMs, put 1 in each of the ***n*** cells. To partition SNP heritability by **functional annotations**, create a SNP info file in which each column (or GRM) represents a category of a functional annotation. Functional annotations may overlap; for example, if a SNP belongs to **6** functional annotation categories, put **1** in each of those **6** cells.

---

## Building GRMs
```sh
for chr in {1..5}
do
  mph --compute_grm --binary_genotype geno --min_maf 0 --snp_info chr.snp_info.csv --snp_weight $chr --num_threads 10 --out $chr
done
```
Create a GRM list, like the example one, **chr.grms.txt**.

---

## Combining GRMs into one
```
mph --merge_grms --grm_list chr.grms.txt --output all_snps
```
If there are two columns in the GRM list file, the second one will be ignored in this procedure. The all-SNPs GRM will be used as the initial value in MINQUE iterations. 

---

## Partitioning SNP heritability
```
mph --minque --binary_grm all_snps --grm_list chr.grms.txt --phenotype phen.csv --trait milk --num_threads 10 --output milk.chr
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
```--constrain``` If enabled, variance component estimates will be constrained to be positive.

```--heritability 0.5```
The SNP heritability value for initializing MINQUE iterations. An accurate value may speed up convergence. The default is 0.5.

```--num_iterations 20```
Max number of MINQUE iterations. The default it 20.

```--rel_tol 1e-4```
Relative tolerence. MINQUE iterations stop when all variance component estimates have a change smaller than that. The default it 1e-4.

```--num_random_vectors 100```
Number of random vectors. The default (100) is usually sufficient. A larger value (such as 500 or 1000) may be needed when convergence issues arise.

---

## Zeroing out small GRM off-diagonal elements
```
mph --zero_grm 0.05 --binary_grm all_snps --output zero_outed
```
All off-diagonal elements smaller than 0.05 are zeroed out, and the resulting matrix is written to the output file.