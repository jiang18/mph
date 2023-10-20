!!! note  
    Utility scripts are hosted in [the MPH repository](https://github.com/jiang18/mph/tree/main/util).

## Utility scripts

### GRM I/O
- `grm_io.R` contains two R functions.
    - write_grm() writes an R GRM into MPH files (.iid and .bin).
    - read_grm() reads an MPH GRM into R.
    - The R matrix to be written can be any **general** relationship matrix.
    - This script is sourced in other R scripts for GRM I/O.
- `grm_txt2bin.R` converts a relationship matrix from txt to MPH binary format. 
    - `Rscript grm_txt2bin.R example.grm.txt example`
    - [example.grm.txt](https://github.com/jiang18/mph/blob/main/examples/example.grm.txt)

### Making files for trait pair


### Heritability enrichment

