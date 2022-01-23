# MINQUE for Partitioning Heritability

## About
MPH is designed to partition SNP heritability with genotypes of related samples or with long-span LDs. For such data, [LDSC](https://github.com/bulik/ldsc), [Haseman-Elston (HE) regression](https://github.com/sriramlab/RHE-mc), or [MQS](https://github.com/genetics-statistics/GEMMA) does not work well.
1. MPH is comparable to GREML in terms of accuracy, while being much faster and more memory-efficient.
2. It can do weighted analyses if residual variances are unequal (such as daughter yield deviations [DYDs] of bulls).
3. It works for many overlapping functional annotations.

---

## Author and contact
[Jicai Jiang](https://cals.ncsu.edu/animal-science/people/jicai-jiang)

---

## Software download
[https://github.com/jiang18/mph/releases/tag/20220120](https://github.com/jiang18/mph/releases/tag/20220120)

## Example data
[The QTL-MAS 2012 data](https://github.com/jiang18/mph/raw/main/QTL-MAS-2012.zip)  
Refer to [https://pubmed.ncbi.nlm.nih.gov/25519515/](https://pubmed.ncbi.nlm.nih.gov/25519515/) for how the data set was simulated.

---

## Input file formats
MPH uses the same file formats as [SSGP](https://sites.google.com/view/ssgp), another program that the same author developed. Refer to [this page](https://sites.google.com/view/ssgp/documentation/manual#h.p_QS3vj5saXQJZ) for details. In short, MPH prefers PLINK bim/fam/bed files for genotypes, and other input files are mostly CSV.

## SNP info file
This CSV file is critical for building genomic relationship matrices (GRMs). Each column corresponds to one GRM and each row corresponds to a SNP. In the example data, there are five columns for five chromosomes.

If a SNP does not belong to a GRM, leave the corresponding cell blank. If a SNP belongs to ***n*** GRMs, put 1 in each of the ***n*** cells. To partition SNP heritability by **functional annotations**, create a SNP info file in which each column (or GRM) represents a category of a functional annotation. Functional annotations may overlap; for example, if a SNP belongs to **6** functional annotation categories, put **1** in each of those **6** cells.

---

## Building GRMs
```sh
for chr in {1..5}
do
  mph --compute_grm --binary_genotype geno --min_maf 0 --min_hwe_pval 1e-8 --snp_info chr.snp_info.csv --snp_weight $chr --num_threads 10 --out $chr
done
```
Create a GRM list, like the example one, **chr.grms.txt**.

---

## Combining GRMs into one
```
mph --merge_grms --grm_list chr.grms.txt --output all_snps
```
If there are two columns in the GRM list file, the second one will be ignored in this procedure. 

---

## Partitioning SNP heritability
```
mph --minque --grm_list chr.grms.txt --phenotype phen.csv --trait milk --num_threads 10 --output milk.chr
```

### To include covariates, add \-\-*covariate_file* and \-\-*covariate_names*.
```
--covariate_file covar.csv --covariate_names all
```
If *all* is specified, MPH will use as covariates all columns from the 2nd to the last in the covariate file.
```
--covariate_file covar.csv --covariate_names g1,g2,g3
```

### For DYD-like data, add \-\-*error_weight_name* to specify individual reliabilies.
```
--error_weight_name milk_wt
```
The error weights can be set to 1/*r*<sup>2</sup>-1.

### Other optional arguments
```--heritability 0.5```
The SNP heritability value for initializing MINQUE iterations. An accurate value may improve convergence. The default is 0.5.

```--num_iterations 20```
Max number of MINQUE iterations. The default it 20.

```--tol 1e-2```
Tolerance. MINQUE iterations stop when logLL has a change smaller than that. The default it 1e-3.

```--num_random_vectors 100```
Number of random vectors. The default (100) is usually sufficient.

---

## Zeroing out small GRM off-diagonal elements
```
mph --zero_grm 0.05 --binary_grm all_snps --output zero_outed
```
All off-diagonal elements smaller than 0.05 are zeroed out, and the resulting matrix is written to the output file.

---

## Simulating phenotypes based on a list of GRMs
```
mph --simulate --num_phenotypes 100 --grm_list chr.grms.txt --heritability 0.5 --output sim_pheno
```
In the *i*th row of the GRM list file are GRM(*i*) and VC(*i*). MPH simulates total genetic values (**g**) by sampling **g** from N(**0**,**V**) in which **V** is equal to the sum of all GRM(*i*)\*VC(*i*). MPH further simulates phenotypes by adding an error term (**e**) to **g** based on heritability.

## Predicting total genetic values
```
mph --pred --mq_file milk.chr --output milk.chr
```
MPH computes total genetic values using the output of --minque and outputs them to a file with a suffix of .mq.gv.csv.
