!!! note  
    Utility scripts are hosted in [the MPH repository](https://github.com/jiang18/mph/tree/main/util).

## Utility scripts

### GRM I/O
MPH computes only **genomic** relationship matrices. `External` **general** relationship matrices can be converted to MPH format for use in MPH. 

- `grm_io.R` contains two R functions.
    - write_grm() writes an R GRM into MPH format files (.iid and .bin).
    - read_grm() reads an MPH GRM into R.
    - The R matrix to be written can be any **general** relationship matrix.
    - This script is sourced in other R scripts for GRM I/O.
- `grm_txt2bin.R` converts a relationship matrix from txt to MPH binary format.
    - `Rscript --no-save grm_txt2bin.R example.grm.txt example`
    - This input file should be structured like [example.grm.txt](https://github.com/jiang18/mph/blob/main/examples/example.grm.txt).
- `AGHmatrix2mph.R` constructs a **numerator** relationship matrix with AGHmatrix and writes it into MPH format files.

### Making files for trait pair
Four R scripts are provided to prepare input files for bivariate REML in MPH.

- 


### Heritability enrichment

