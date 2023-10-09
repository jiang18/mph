## GRM

### Input file formats
MPH uses PLINK bim/fam/bed files for genotypes. Other input files are comma-separated values (CSV) or space-delimited text.

### SNP info file
This CSV file is necessary for making genomic relationship matrices (GRMs). 

The file has a header line. The first column lists SNP IDs. Each of other columns corresponds to one GRM and each row corresponds to a SNP. In the [example data](examples.md#qtl-mas-2012), [**chr.snp_info.csv**](https://github.com/jiang18/mph/blob/main/examples/QTL-MAS-2012/chr.snp_info.csv?plain=1) has five columns for five chromosomes.

If a SNP does not belong to a GRM, leave the corresponding cell **empty**. If a SNP belongs to ***n*** GRMs, put **1** in each of the ***n*** cells. To partition SNP heritability by **functional annotations**, create a SNP info file in which each column (or GRM) represents a category of a functional annotation. Functional annotations may overlap; for example, if a SNP belongs to **6** functional annotation categories, put **1** in each of those **6** cells.

### Making additive GRMs
```sh
for chr in {1..5}
do
  mph --make_grm --binary_genotype geno --min_maf 0 --min_hwe_pval 1e-8 --snp_info chr.snp_info.csv --snp_weight $chr --num_threads 10 --out $chr
done
```

### Making dominance GRMs
To construct a dominance GRM, add `--dominance`.
```sh
for chr in {1..5}
do
  mph --make_grm --dominance --binary_genotype geno --min_maf 0 --min_hwe_pval 1e-8 --snp_info chr.snp_info.csv --snp_weight $chr --num_threads 10 --out $chr.dom
done
```

### GRM list
The GRM list file is a space-delimited text file listing GRMs, like the example one, [**chr.grms.txt**](https://github.com/jiang18/mph/blob/main/examples/QTL-MAS-2012/chr.grms.txt). It has no header and one or two columns. The first column must be the file path for each GRM. The second column is optional and can be a number or a short GRM identifier. 

The GRM list is needed for `--merge_grms`, `--deduct_grms`, `--make_core`, `--make_fore`, and `--minque`.

### Merging GRMs
If a SNP is used in multiple GRMs, the SNP will be treated to be multiple identical SNPs in the resulting GRM.
```sh
mph --merge_grms --grm_list chr.grms.txt --output all_snps
```
If there are two columns in the GRM list file, the second one will be ignored in this procedure. 

### Deducting GRMs
```sh
mph --deduct_grms --grm_list list.grms.txt --output deducted
```

### Zeroing out GRM elements
```sh
mph --zero_grm 0.05 --binary_grm all_snps --output zero_outed
```
All off-diagonal elements smaller than 0.05 are zeroed out, and the resulting matrix is written to the output file.

## REML/MINQUE

### Input
- `--grm_list`: a space-delimited text file without header.
    - The first column lists GRM file path. 
    - The second column is optional and can list an initial VC value for each GRM.
- `--phenotype`: a CSV file with a header line.
    - The first column must be the individual ID.
- `--covariate_file`: a CSV file with a header line.
    - The first column must be the individual ID.
    - The covariate file is optional.

Missing values of phenotypes or covariates need to be left empty. Do not use space, -9, NA, or NaN.

### REML or iterative MINQUE
```
--minque
```
This flag turns on iterative MINQUE (equivalent to the Fisher scoring algorithm for REML).

```sh
mph --minque --grm_list chr.grms.txt --phenotype phen.csv --trait milk --num_threads 10 --output milk.chr
```

### MINQUE(0) or MINQUE(1)
To force MINQUE(0) or MINQUE(1), add `--num_iterations 1`. The second column of the GRM list file should be set to 0 for MINQUE(0) and 1 for MINQUE(1).
```sh
mph --minque --num_iterations 1 --grm_list chr.grms.txt --phenotype phen.csv --trait milk --num_threads 10 --output milk.chr
```

### Memory-saving mode
```
--save_memory
```
This flag enables the memory-saving mode. This mode is not necessarily slower, particularly on a solid-state drive.

### Covariates
To include covariates, add `--covariate_file` and `--covariate_names`.
```sh
--covariate_file covar.csv --covariate_names all
```
If **all** is specified, MPH will use as covariates all columns from the 2nd to the last in the covariate file.
```sh
--covariate_file covar.csv --covariate_names g1,g2,g3
```

### Pseudo-phenotypes
If needed, add `--error_weight_name` to specify individual reliabilies for pseudo-phenotypes (e.g., de-regressed estimated breeding values).
```
--error_weight_name milk_wt
```
The error weights can be set to 1/*r*<sup>2</sup>-1.

### General relationship matrix
The `--grm_list` file can include any **general** relationship matrix, not just **genomic** relationship matrices.

### Additional options
```--heritability 0.5```
The heritability value for initializing REML/MINQUE iterations. An accurate value may improve convergence. The default is 0.5.

```--num_iterations 20```
Max number of REML iterations. The default is 20.

```--tol 0.01```
Tolerance. REML iterations stop when logLL has a change smaller than that. The default is 0.01.

```--num_random_vectors 100```
Number of random probing vectors for stochastic trace estimation. The default (100) is usually sufficient.

```--seed 0```
Random seed in stochastic trace estimation. The default is 0.

### Output
| Filename suffix | Description |
|----------|----------|
| mq.blue.csv | Best linear unbiased estimates of fixed effects |
| mq.iter.csv | Summary of REML iterations |
| mq.vc.csv | Variance component estimates |
| mq.py.csv | Residuals of the linear mixed model |

**mq.vc.csv** has the following columns.

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

## CORE/FORE

### COvariance between Random Effects (CORE)

### First-Order interaction between Random Effects (FORE)

## Simulation
Simulating phenotypes based on a list of GRMs
```sh
mph --simulate --num_phenotypes 100 --grm_list chr.grms.txt --heritability 0.5 --output sim_pheno
```
In the *i*th row of the GRM list file are GRM(*i*) and VC(*i*). MPH simulates total genetic values (**g**) by sampling **g** from N(**0**,**V**) in which **V** is equal to the sum of all GRM(*i*)\*VC(*i*). MPH further simulates phenotypes by adding an error term (**e**) to **g** based on heritability.

## Prediction
Empirical best linear unbiased predictions (EBLUPs)
```sh
mph --pred --mq_file milk.chr --output milk.chr
```
MPH computes EBLUPs using the output of `--minque` and outputs them to a file with a suffix of .mq.blup.csv. For genomic partitioning, EBLUPs are the estimates of direct genomic values. 

