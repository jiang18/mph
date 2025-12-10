!!! note  
    Option names may be abbreviated if the abbreviation is unique or is an exact match for some defined option; e.g., `--phenotype` works the same as `--phenotype_file`.

## Making a GRM from SNPs

### Input
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--bfile`<br>`--binary_genotype_file` | FILE&nbsp;PREFIX | Required | PLINK bed/bim/fam filename prefix |
| `--snp_info_file` | FILE | Required | SNP info file in CSV format |
| `--snp_weight_name` |  STRING | Optional | Specify a column name in the SNP info file for weighting SNPs in a GRM. If not set, all variants in the SNP info file will be used with a default weight of 1. |
| `--min_maf` | FLOAT | Optional | Filter out all variants with minor allele frequency (MAF) less than or equal to the provided threshold [default=0] |
| `--min_hwe_pval` | FLOAT | Optional | Filter out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold [default=0] |

#### Genotype files
MPH uses PLINK bim/fam/bed files for genotypes. The bed/bim/fam file format is described on the [PLINK website](https://www.cog-genomics.org/plink/1.9/formats).

#### SNP info file
The SNP info file is a CSV file with a header line. The first column must list SNP IDs. Each of the other columns corresponds to one GRM, and each row corresponds to a SNP.

If a SNP does not belong to a GRM, leave the corresponding cell **empty**. If a SNP belongs to ***n*** GRMs, put **1** in each of the ***n*** cells. To partition SNP heritability by **functional annotations**, create a SNP info file in which each column (or GRM) represents a functional annotation category. Functional annotations may overlap; for example, if a SNP belongs to **6** functional annotation categories, put **1** in each of those **6** cells.

In the example of [partitioning heritability by chromosomes](examples.md#by-chromosomes), [**chr.snp_info.csv**](https://github.com/jiang18/mph/blob/main/examples/QTL-MAS-2012/chr.snp_info.csv?plain=1) has five columns for five chromosomes.

[A Perl script](scripts.md#making-a-snp-info-file) is provided to create the SNP info file for partitioning heritability across functional annotations.

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

The GRM calculation generates two files:

1. A file with the suffix `.grm.iid`
2. A file with the suffix `.grm.bin`

For dominance GRMs, the suffixes are modified as follows:

1. `.d.grm.iid`
2. `.d.grm.bin`

### Custom genotype coding
MPH includes a feature for computing GRM with custom genotype coding.

To use custom genotype coding:

1. Add 3 columns to the [SNP info file](#snp-info-file) specifying the numeric codes for A<sub>1</sub>A<sub>1</sub>, A<sub>1</sub>A<sub>2</sub>, and A<sub>2</sub>A<sub>2</sub> genotypes for each genomic variant.
2. Name these columns (e.g., header1, header2, header3) in the header line.
3. When running `mph --make_grm`, include the option `--snp_genotype_coding header1,header2,header3` to enable GRM computation using the custom genotype coding specified in the SNP info file.

Important notes:

- The `--snp_genotype_coding` option requires exactly 3 tokens. Ensure that the first, second, and third token match A<sub>1</sub>A<sub>1</sub>, A<sub>1</sub>A<sub>2</sub>, and A<sub>2</sub>A<sub>2</sub>, respectively.
- A<sub>1</sub> and A<sub>2</sub> refer to the alleles in columns 5 and 6 of the [PLINK .bim file](https://www.cog-genomics.org/plink/1.9/formats#bim), respectively.
- To provide users with complete control over GRM computation, custom genotype codes are used as provided, without any centering applied.
- The `--snp_weight_name` option remains available for weighting genotype codes.

## Making a GRM from GRMs

### Input
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--grm_list` | FILE | Required | Space-delimited text file listing GRMs (include paths if needed) |

The GRM list file is a space-delimited text file, such as [**AD.grms.txt**](https://github.com/jiang18/mph/blob/main/examples/QTL-MAS-2012/AD.grms.txt). It has no header and one or two columns. 
The first column must be the file path for each GRM. The second column is optional and can be a short GRM identifier for `--make_fore` and `--make_core`. 

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

## Subsetting a GRM
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--subset_grm` | FLAG | Required | Flag to initialize the subsetting routine |
| `--binary_grm_file` | FILE PREFIX | Required | Filename prefix of input GRM to be subsetted |
| `--keep` | FILE | Required | Text file listing individual IDs to keep in the subset |
| `--output_file` | FILE PREFIX | Required | Output filename prefix |

The file specified by `--keep` is a text file with one individual ID per line and no header. It determines both which individuals to keep and their order in the output GRM.

The program stops with an error message if any individual IDs listed in the `--keep` file are not found in the input GRM. The output GRM will have individuals in the exact order specified in the `--keep` file.

## REML/MINQUE

### Input
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--grm_list`<br>`--binary_grm_file` | FILE<br>FILE&nbsp;PREFIX | Required | Use one or both: space-delimited text file listing GRMs (`--grm_list`) and/or filename prefix of a single binary GRM (`--binary_grm_file`); include paths if needed |
| `--phenotype_file` | FILE | Required | Phenotype file in CSV format |
| `--trait` | STRING | Required | Single trait name or comma-separated list of trait names that should match the column names in the phenotype file |
| `--error_weight_names` | STRING | Optional | Single column name or comma-separated list of column names in the phenotype file, specifying individual error variance weights for each corresponding trait |
| `--covariate_file` | FILE | Optional | Covariate file in CSV format |
| `--covariate_names` | STRING | Optional | Comma-separated list of covariates to include in the analysis |

#### GRM input
Specify GRMs using one or both of the following options:

- `--grm_list` with a list file containing one or more GRMs
- `--binary_grm_file` to directly specify a single GRM

If both options are provided, the GRM specified by `--binary_grm_file` will be appended to the list from `--grm_list`.

The GRM list file is a space-delimited text file, such as [**chr.grms.txt**](https://github.com/jiang18/mph/blob/main/examples/QTL-MAS-2012/chr.grms.txt). It has no header and one or two columns. 
The first column must be the file path for each GRM. The second column is optional (defaults to 1) and can list an initial variance component (VC) value for each GRM. 
When `--binary_grm_file` is used with `--grm_list`, its initial VC value is set to the average of the VC values from the list file.

#### Phenotype file
The phenotype file is a CSV file with a header line. The first column must be the individual ID. The header line needs to contain trait names.

Missing values of phenotypes need to be left empty. **Do not use space, -9, NA, or NaN.**

MPH has `--error_weight_names` to accommodate individual reliabilities (*r*<sup>2</sup>) for pseudo-phenotypes (e.g., de-regressed estimated breeding values). The error weights can be set to 1/*r*<sup>2</sup>-1 and kept as a column of the phenotype file.

#### Covariate file
The covariate file is a CSV file with a header line. The first column must be the individual ID. The header line needs to contain covariate names. All covariate values must be numeric.

Missing values of covariates need to be left empty. **Do not use space, -9, NA, or NaN.**

If `--covariate_names all` is specified, MPH will use all columns (except the first) as covariates.

*Intercept handling:*
- If `--covariate_names` is not specified, MPH automatically includes an intercept term.
- If `--covariate_names` is specified, MPH does not automatically add an intercept. In this case, users must include a column of 1's in the covariate file and specify it in `--covariate_names` if an intercept is desired.

### Options
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--minque`<br>`--reml` | FLAG | Required | Flag to run REML or iterative MINQUE. This option can be specified as `--minque` or `--reml`, and they are interchangeable. |
| `--save_memory` | FLAG | Optional | Flag to enable the memory-saving mode |
| `--num_threads` | INT | Optional | Number of computational threads to use [default=1] |
| `--heritability` | FLOAT | Optional | Initial heritability value(s) in REML iterations [default=0.5] |
| `--num_grms` | INT | Optional | Number of the list's first matrices to be used as genomic relationships [default=all] |
| `--num_iterations` | INT | Optional | Max number of REML iterations [default=20] |
| `--tol` | FLOAT | Optional | Absolute convergence tolerance for REML log-likelihood [default=0.01] |
| `--num_random_vectors` | INT | Optional | Number of random probing vectors for stochastic trace estimation [default=100] |
| `--distribution` | STRING | Optional | Specify whether the stochastic trace estimation uses a **Gaussian** or **Rademacher** distribution [default=Rademacher] |
| `--seed` | INT | Optional | Random seed in stochastic trace estimation [default=0] |

MPH incorporates Fisher's scoring for REML, equivalent to iterative MINQUE, alongside a trust-region dogleg algorithm. `--minque` mirrors `--reml`.

To force MINQUE(0) or MINQUE(1), set `--num_iterations 1`. The second column of the GRM list file should be set to 0 for MINQUE(0) and 1 for MINQUE(1).

The memory-saving mode (`--save_memory`) is not necessarily slower, especially when file caching is effectively utilized. 

### Output
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--output_file` | FILE PREFIX | Required | Output filename prefix |

`mph --reml` produces four or five output files, as shown below.

| Filename suffix | Description |
|----------|----------|
| .mq.blue.csv | Best linear unbiased estimates of fixed effects |
| .mq.iter.csv | Summary of REML iterations |
| .mq.vc.csv | Variance component estimates |
| .mq.py.csv | Residuals of the linear mixed model |
| .mq.cor.csv | Estimates and standard errors of correlations between traits. **Exclusively available for multi-trait analyses.** |

**.mq.vc.csv** has the following columns.

| Column | Description |
|----------|----------|
| trait_x | Trait name |
| trait_y | Trait name. If the same as trait_x, components correspond to the trait variance; otherwise, the covariance between trait_x and trait_y. |
| vc_name | Name of the (co)variance component or GRM |
| m | Weighted number of SNPs in the GRM |
| var | Estimate of the (co)variance component |
| seV | Standard error of the (co)variance component estimate |
| pve | Estimate of the proportion of (co)variance explained |
| seP | Standard error of the estimate of the proportion of (co)variance explained |
| enrichment | Estimate of the per-SNP heritability enrichment |
| seE | Standard error of the estimate of the per-SNP heritability enrichment |

Additional columns display sampling covariance matrices of estimates for enrichments and for all (co)variance components.

!!! note  
    Estimates of PVEs and enrichments are valid only when functional annotation categories do not overlap with one another. If functional categories actually overlap, one more quick computation is needed to [recompute PVEs and enrichments](scripts.md#from-vcs-to-enrichments). 

## Simulation
Simulating phenotypes based on a list of GRMs
```sh
mph --simulate --num_phenotypes 100 --grm_list chr.grms.txt --heritability 0.5 --output pheno
```
In the *i*th row of the GRM list file are GRM(*i*) and VC(*i*). MPH simulates total genetic values (**g**) by sampling **g** from N(**0**,**V**) in which **V** is equal to the sum of all GRM(*i*)\*VC(*i*). MPH further simulates phenotypes by adding an error term (**e**) to **g** based on heritability.

## Prediction
Empirical best linear unbiased predictions (EBLUPs)
```sh
mph --pred --mq_file milk --output milk
```
MPH computes EBLUPs using the output of `--reml` and outputs them to a file with a suffix of **.mq.blup.csv**. For genomic partitioning, EBLUPs are the estimates of direct genomic values. 

**Multi-trait BLUP is not currently supported.**

## General relationship matrix
The `--grm_list` file can list any **general** relationship matrix, not necessarily a **genomic** relationship matrix.

[Auxiliary scripts](scripts.md#grm-inputoutput) are provided to convert a general relationship matrix to the MPH format.
