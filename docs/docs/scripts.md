!!! note  
    Auxiliary scripts are hosted in [the MPH repository](https://github.com/jiang18/mph/tree/main/scripts).

## GRM input/output
MPH computes only **genomic** relationship matrices. **External general** relationship matrices can be converted to MPH format for use in MPH. 

- [**mph_functs.R**](https://github.com/jiang18/mph/blob/main/scripts/mph_functs.R) contains two R functions for GRM input/output.
    - `write_grm()` writes an R GRM into MPH format files (.iid and .bin).
    - `read_grm()` reads an MPH GRM into R.
    - The R matrix to be written can be any **general** relationship matrix.
    - This script is sourced in other R scripts for GRM input/output.
- [**grm_txt2bin.R**](https://github.com/jiang18/mph/blob/main/scripts/grm_txt2bin.R) converts a relationship matrix from txt to MPH binary format.
    - `Rscript --no-save grm_txt2bin.R example.grm.txt example`
    - This input file should be structured like [**example.grm.txt**](https://github.com/jiang18/mph/blob/main/examples/example.grm.txt).
- [**AGHmatrix2mph.R**](https://github.com/jiang18/mph/blob/main/scripts/AGHmatrix2mph.R) constructs a **numerator** relationship matrix with [AGHmatrix](https://cran.r-project.org/web/packages/AGHmatrix/) and writes it into MPH format files.
- [**make_gci_grm.R**](https://github.com/jiang18/mph/blob/main/scripts/make_gci_grm.R) generates a genotype–covariate interaction GRM from a routine GRM and an indicator matrix.

## Heritability enrichment
### Making a SNP info file
A Perl script, [**make_snp_info.pl**](https://github.com/jiang18/mph/blob/main/scripts/make_snp_info.pl), is provided to create [the SNP info file](options.md#snp-info-file) for partitioning heritability across functional annotations.

Usage:

- `perl make_snp_info.pl PLINK-bim-file functional-annotation-file output-filename-prefix`
- Both input files should be in plain text format without header lines.
- Functional-annotation-file's first four columns should be chrom, start, end, and category.

Example:

- `perl make_snp_info.pl test.bim test.annot.txt test`
- The input files are [**test.bim**](https://github.com/jiang18/mph/blob/main/examples/test.bim) and [**test.annot.txt**](https://github.com/jiang18/mph/blob/main/examples/test.annot.txt).
- An output file named **test.snp_info.csv** will be generated. 

### From VCs to enrichments
The **.mq.vc.csv** file produced by `mph --reml` contains the estimates and SEs of VCs, PVEs, and enrichments, where PVEs and enrichments are valid only when functional annotation categories do not overlap with one another.

If functional categories actually overlap, one more quick computation is needed to recompute PVEs and enrichments from **.mq.vc.csv** and [the SNP info file](options.md#snp-info-file). This can be quickly done by `recompute_enrichments()` in [**mph_functs.R**](https://github.com/jiang18/mph/tree/main/scripts/mph_functs.R).

`recompute_enrichments(vcfile, crossprod, nsnps=NA, annot.size=NA, trait.x=NA, trait.y=NA)`

- vcfile: the .mq.vc.csv file produced by `mph --reml`.
- crossprod: the crossproduct of SNP incidence matrix and SNP weighting matrix. Its row names and columns should match the annotation categories of interest and the rows of vcfile, respectively.
- nsnps: the total number of SNPs. If `NA`, it is set to the first GRM's number of SNPs in vcfile.
- annot.size: a list of the number of SNPs in each annotation category listed in the row names of crossprod. It can be computed from the column-wise sum of the corresponding SNP incidence matrix. If `NA`, it is set to the `m` column of vcfile.
- trait.x and trait.y: trait names matching the first two columns of vcfile. If both are `NA`, they are set to the trait names in the first line of vcfile. Otherwise, if only one is `NA`, it is set to the value of the other.

A usage example is available [here](examples.md#by-functional-annotations).
