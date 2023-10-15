## Simulated datasets

### QTL-MAS 2012
- The data set is available for download [here](https://github.com/jiang18/mph/raw/main/examples/QTL-MAS-2012.zip).
- [This article](https://bmcproc.biomedcentral.com/articles/10.1186/1753-6561-8-S5-S1) describes how the data set was simulated.
- 3k related animals, 10k SNPs, 5 chromosomes, and 3 traits

### Sequence genotypes
- The data set is available for download [here]().
- Sequence genotypes were simulated by [genosim](https://aipl.arsusda.gov/software/genosim/).
- 10k unrelated individuals, 5M sequence variants, and 30 chromosomes
- Functional annotations were quickly simulated by assigning [LDSC baseline annotations](https://console.cloud.google.com/storage/browser/_details/broad-alkesgroup-public-requester-pays/LDSCORE/1000G_Phase3_baseline_ldscores.tgz) to the simulated sequence variants in order.
- Phenotypes were simulated using the [S-LDSC baseline model enrichment estimates for human traits](https://www.nature.com/articles/ng.3404/figures/4).
    - The enrichment estimates were used to compute variance component (VC) estimates for intercept and 24 main functional annotations.
    - The VC estimates were set as true values in [phenotype simulation](options.md#simulation).
    - A small value was added to intercept's VC to enhance the all-in-one GRM's positive definiteness. 

## Partitioning heritability

### By chromosomes
Partitioning heritability by chromosomes for the [QTL-MAS 2012](#qtl-mas-2012) dataset

1. Create a [SNP info file](options.md#snp-info-file): [**chr.snp_info.csv**](https://github.com/jiang18/mph/blob/main/examples/QTL-MAS-2012/chr.snp_info.csv?plain=1).
2. [Make a GRM](options.md#making-a-grm-from-snps) for each chromosome
    ```sh
    mkdir grm

    for chr in {1..5}
    do
    mph --make_grm --binary_genotype geno --min_maf 0 --min_hwe_pval 1e-8 --snp_info chr.snp_info.csv --snp_weight $chr --num_threads 10 --out ./grm/$chr
    done
    ```
3. Create a [GRM list](options.md#grm-list-file): [**chr.grms.txt**](https://github.com/jiang18/mph/blob/main/examples/QTL-MAS-2012/chr.grms.txt).
4. Run [REML/MINQUE](options.md#remlminque).
```sh
mkdir minque

mph --minque --grm_list chr.grms.txt --phenotype phen.csv --trait milk --error_weight milk_wt --covariate_file covar.csv --covariate_names all --num_threads 10 --out ./minque/milk.chr
```

### By functional annotations

## Dominance and epistasis

## Genetic correlation
