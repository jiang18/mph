## Make additive GRM
- PLINK bed/bim/fam: geno
- SNP info file: chr.snp_info.csv

```sh
mkdir grm

for chr in {1..5}
do
  mph --make_grm --binary_genotype geno --min_maf 0 --min_hwe_pval 1e-8 --snp_info chr.snp_info.csv --snp_weight $chr --num_threads 10 --out ./grm/$chr
done
```

## REML (or I-MINQUE)
- GRM list file: chr.grms.txt
- Phenotype file: phen.csv
- Covariate file: covar.csv

```sh
mkdir minque

mph --minque --grm_list chr.grms.txt --phenotype phen.csv --trait milk --error_weight milk_wt --covariate_file covar.csv --covariate_names all --num_threads 10 --out ./minque/milk.chr
```
