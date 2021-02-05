# MPH: MINQUE for Partitioning Heritability
## Software download
https://github.com/jiang18/mph/releases/tag/20210204

## Example: the QTL-MAS 2012 data
Refer to https://bmcproc.biomedcentral.com/articles/10.1186/1753-6561-8-S5-S1 for how the data set was simulated.

Refer to https://sites.google.com/view/ssgp/documentation/manual#h.p_QS3vj5saXQJZ for the data format.

MPH prefers Plink bim/fam/bed files for genotypes. Other input files are mostly CSV.

## SNP info file
This CSV file is critical for building genomic relationship matrices (GRM). Each column corresponds to one GRM. In the example data, there are five columns for five chromosomes.

If a SNP does not belong to a GRM, leave the corresponding cell blank. If a SNP belongs to ***n*** GRMs, put **1/*n*** in the ***n*** cells. The row sum of any SNP should be equal to 1.

To partition SNP heritability by **functional annotations**, create a SNP info file in which each column (or GRM) represents a category of a functional annotation. Functional annotations may overlap. If a SNP belongs to ***n*** functional annotation categories, put **1/*n*** in those ***n*** cells.

## Building GRMs
```sh
for chr in {1..5}
do
  mph --compute_grm --binary_genotype geno --min_maf 0 --snp_info chr.snp_info.csv --snp_weight $chr --num_threads 10 --out $chr
done
```
Create a GRM list, like the example one, **chr.grms.txt**.

## Combining GRMs into one
```
mph --combine_grms --grm_list chr.grms.txt --output allSnps
```
The all-SNPs GRM will be used as the initial value in MINQUE iterations.

## Partitioning SNP heritability
```
mph --minque --binary_grm allSNPs --grm_list chr.grms.txt --phenotype phen.csv --trait milk --num_threads 10 --output milk.chr
```

### To include covariates, add *--covariate_file* and *--covariate_names*.
```
--covariate_file covar.csv --covariate_names all
```
```
--covariate_file covar.csv --covariate_names g1,g2,g3
```

### For DYD or similar data, MPH has an argument (*--error_weight_name*) to specify individual accuracies.
```
--error_weight_name milk_wt
```
The error weights can be set to 1/r^2-1.

### Other optional arguments
```--heritability 0.5```
The SNP heritability value for initializing MINQUE iterations. An accurate value may speed up convergence. The default is 0.5.

```--num_iterations 20```
Max number of MINQUE iterations. The default it 20.

```--rel_tol 1e-5```
Relative tolerence. MINQUE iterations stop when all variance component estimates have a change smaller than that. The default it 1e-5.

```--num_rademacher 100```
Number of Rademacher samples. The default (100) is generally sufficient. A larger value (500 or 1000) can work better for small samples (<1000).

