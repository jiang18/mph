---
title: GRM
layout: template
filename: grm.md
---

## Input file formats
MPH uses the same file formats as [SSGP](https://sites.google.com/view/ssgp), another program that the same author developed. Refer to [this page](https://sites.google.com/view/ssgp/documentation/manual#h.p_QS3vj5saXQJZ) for details. In short, MPH prefers PLINK bim/fam/bed files for genotypes, and other input files are mostly CSV.

## SNP info file
This CSV file is critical for making genomic relationship matrices (GRMs). Each column corresponds to one GRM and each row corresponds to a SNP. In the example data, there are five columns for five chromosomes.

If a SNP does not belong to a GRM, leave the corresponding cell blank. If a SNP belongs to ***n*** GRMs, put 1 in each of the ***n*** cells. To partition SNP heritability by **functional annotations**, create a SNP info file in which each column (or GRM) represents a category of a functional annotation. Functional annotations may overlap; for example, if a SNP belongs to **6** functional annotation categories, put **1** in each of those **6** cells.

---

## Making additive GRMs
```sh
for chr in {1..5}
do
  mph --make_grm --binary_genotype geno --min_maf 0 --min_hwe_pval 1e-8 --snp_info chr.snp_info.csv --snp_weight $chr --num_threads 10 --out $chr
done
```

## Making dominance GRMs
To construct a dominance GRM, add --dominance.
```sh
for chr in {1..5}
do
  mph --make_grm --dominance --binary_genotype geno --min_maf 0 --min_hwe_pval 1e-8 --snp_info chr.snp_info.csv --snp_weight $chr --num_threads 10 --out $chr.dom
done
```

---

## Combining GRMs into one
If a SNP is used in multiple GRMs, the SNP will be treated to be multiple identical SNPs in the resulting GRM.
```
mph --merge_grms --grm_list chr.grms.txt --output all_snps
```
If there are two columns in the GRM list file, the second one will be ignored in this procedure. 

## Deducting GRMs

## Zeroing out small GRM off-diagonal elements
```
mph --zero_grm 0.05 --binary_grm all_snps --output zero_outed
```
All off-diagonal elements smaller than 0.05 are zeroed out, and the resulting matrix is written to the output file.

---

## GRM list
Create a GRM list, like the example one, **chr.grms.txt**.
