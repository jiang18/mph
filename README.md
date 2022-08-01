# MPH
## MINQUE for Partitioning Heritability

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
[https://github.com/jiang18/mph/releases/tag/20220729](https://github.com/jiang18/mph/releases/tag/20220729)

## Example data
[The QTL-MAS 2012 data](https://github.com/jiang18/mph/raw/main/QTL-MAS-2012.zip)  
Refer to [https://pubmed.ncbi.nlm.nih.gov/25519515/](https://pubmed.ncbi.nlm.nih.gov/25519515/) for how the data set was simulated.

---

## [Genomic relationship matrix](./grm.md)

---

## [Genomic partitioning (MINQUE/REML)](./minque.md)

---

## [Simulation of phenotypes](./simulation.md)

---

## [Empirical best linear unbiased predictions](./prediction.md)
