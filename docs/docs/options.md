!!! note
    Option names may be abbreviated if the abbreviation is unique or is an exact match for some defined option; e.g., `--phenotype` works the same as `--phenotype_file`.

## Making a GRM from genotypes

### Input
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--binary_genotype_file` | FILE PREFIX | Required | PLINK bed/bim/fam filename prefix |
| `--snp_info_file` | FILE | Required | SNP info file in CSV format |
| `--snp_weight_name` |  STRING | Optional | Must be a column header in the SNP info file and specify which column to use for weighting SNPs in a GRM. If not set, all variants in the SNP info file will be used with a weight of 1. |
| `--min_maf` | FLOAT | Optional | Filter out all variants with minor allele frequency (MAF) less than or equal to the provided threshold [default=0] |
| `--min_hwe_pval` | FLOAT | Optional | Filter out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold [default=0] |

#### Genotype files
MPH uses PLINK bim/fam/bed files for genotypes. The bed/bim/fam file format is described on the [PLINK website](https://www.cog-genomics.org/plink/1.9/formats).

#### SNP info file
This CSV file must have a header line. The first column lists SNP IDs. Each of other columns corresponds to one GRM and each row corresponds to a SNP. In the example of [partitioning heritability by chromosomes](examples.md#qtl-mas-2012), [**chr.snp_info.csv**](https://github.com/jiang18/mph/blob/main/examples/QTL-MAS-2012/chr.snp_info.csv?plain=1) has five columns for five chromosomes.

If a SNP does not belong to a GRM, leave the corresponding cell **empty**. If a SNP belongs to ***n*** GRMs, put **1** in each of the ***n*** cells. To partition SNP heritability by **functional annotations**, create a SNP info file in which each column (or GRM) represents a category of a functional annotation. Functional annotations may overlap; for example, if a SNP belongs to **6** functional annotation categories, put **1** in each of those **6** cells.

### Options
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--make_grm` | FLAG | Required | To make a GRM from genotypes |
| `--dominance` | FLAG | Optional | Flag to make a dominance GRM rather than an additive GRM |
| `--num_threads` | INT | Optional | Number of computational threads to use [default=1] |

### Output
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--output_file` | FILE PREFIX | Required | Output filename prefix |

Two files will be generated for a GRM: one with the suffix **.grm.iid** and the other with **.grm.bin**.

## Making a GRM from GRMs

### Input
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--grm_list` | FILE | Required | Space-delimited text file listing GRMs that have been made |

The GRM list file is a space-delimited text file, such as [**AD.grms.txt**](https://github.com/jiang18/mph/blob/main/examples/QTL-MAS-2012/AD.grms.txt). It has no header and one or two columns. The first column must be the file path for each GRM. The second column is optional and can be a short GRM identifier for `--make_fore` and `--make_core`. 

### Options
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--merge_grms` | FLAG | Optional | Flag to merge a list of GRMs |
| `--deduct_grms` | FLAG | Optional | Flag to deduct the subsequent GRMs from the first in a list |
| `--make_fore` | FLAG | Optional | Flag to make first-order interaction GRMs |
| `--make_core` | FLAG | Optional | Flag to make GRMs for covariances between random effects |
| `--num_threads` | INT | Optional | Number of computational threads to use [default=1] |

If a SNP is used in multiple GRMs, the SNP will be treated to be multiple identical SNPs in the resulting GRM of `--merge_grms`.

### Output
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--output_file` | FILE PREFIX | Required | Output filename prefix |

- `--merge_grms` and `--deduct_grms` generate one GRM.
- `--make_fore` generates *n*(*n*+1)/2 GRMs for *n* GRMs in a GRM list.
- `--make_core` generates *n*(*n*-1)/2 GRMs for *n* GRMs in a GRM list.

## REML/MINQUE

### Input
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--grm_list` | FILE | Required | Space-delimited text file listing GRMs that have been made |
| `--phenotype_file` | FILE | Required | Phenotype file in CSV format |
| `--trait` | STRING | Required | Must be a column header in the phenotype file and specify which trait to analyze |
| `--error_weight_name` | STRING | Optional | Must be a column header in the phenotype file and specify which column contains individual error variance weights |
| `--covariate_file` | FILE | Optional | Covariate file in CSV format |
| `--covariate_names` | STRING | Optional | Comma-separated list of covariates to include in the analysis |

#### GRM list file
The GRM list file is a space-delimited text file, such as [**chr.grms.txt**](https://github.com/jiang18/mph/blob/main/examples/QTL-MAS-2012/chr.grms.txt). It has no header and one or two columns. The first column must be the file path for each GRM. The second column is optional and can list an initial variance component (VC) value for each GRM. 

#### Phenotype file
The phenotype file is a CSV file with a header line. The first column must be the individual ID. The header line needs to contain trait names.

Missing values of phenotypes need to be left empty. Do not use space, -9, NA, or NaN.

MPH has `--error_weight_name` to accommodate individual reliabilies (*r*<sup>2</sup>) for pseudo-phenotypes (e.g., de-regressed estimated breeding values). The error weights can be set to 1/*r*<sup>2</sup>-1 and kept as a column of the phenotype file.

#### Covariate file
The covariate file is a CSV file with a header line. The first column must be the individual ID. The header line needs to contain covariate names. 

Missing values of covariates need to be left empty. Do not use space, -9, NA, or NaN.

If `--covariate_names all` is specified, MPH will use as covariates all columns from the 2nd to the last in the covariate file.

### Options
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--minque` | FLAG | Required | Flag to run REML or iterative MINQUE |
| `--save_memory` | FLAG | Optional | Flag to enable the memory-saving mode |
| `--num_threads` | INT | Optional | Number of computational threads to use [default=1] |
| `--heritability` | FLOAT | Optional | Initial heritability value in REML iterations [default=0.5] |
| `--exact_init` | FLAG | Optional | Flag to set initial VC values to the exact values of the GRM list file's second column rather than their ratios |
| `--num_iterations` | INT | Optional | Max number of REML iterations [default=20] |
| `--tol` | FLOAT | Optional | Absolute convergence tolerance for REML log-likelihood [default=0.01] |
| `--num_random_vectors` | INT | Optional | Number of random probing vectors for stochastic trace estimation [default=100] |
| `--seed` | INT | Optional | Random seed in stochastic trace estimation [default=0] |

To force MINQUE(0) or MINQUE(1), set `--num_iterations 1`. The second column of the GRM list file should be set to 0 for MINQUE(0) and 1 for MINQUE(1).

The memory-saving mode (--save_memory) is not necessarily slower, particularly on a solid-state drive. 

### Output
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--output_file` | FILE PREFIX | Required | Output filename prefix |

Four files are generated, as shown below.

| Filename suffix | Description |
|----------|----------|
| .mq.blue.csv | Best linear unbiased estimates of fixed effects |
| .mq.iter.csv | Summary of REML iterations |
| .mq.vc.csv | Variance component estimates |
| .mq.py.csv | Residuals of the linear mixed model |

**.mq.vc.csv** has the following columns.

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

## General relationship matrix
The `--grm_list` file can list any **general** relationship matrix, not necessarily a **genomic** relationship matrix.

[Scripts](util.md) are provided to convert a general relationship matrix to the MPH format. 

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
