# Clean the Familinx pedigree
```sh
mkdir familinx
cd familinx
perl clean_familinx_step1.pl
perl clean_familinx_step2.pl
cd ..
```
# Extract a subset from the Familinx pedigree
```sh
perl subset_ped.pl 1
```
# Simulate genotypes based on the extracted pedigree
```sh
cd 1
../markersim
../genosim
perl ../aipl2plink.pl 10k
plink --file 10k --make-bed --out 10k
```
# GRM
```sh
plink --bfile 10k --make-rel square --out grm
```
```{r}
# Check relatedness
library(data.table)
G = fread("grm.rel")
G = as.matrix(G)
summary(G[lower.tri(G)])
```
# LD r2
```sh
plink --bfile 10k --r2 square --out r2
```
```{r}
# Check inter-chromosome LD r2
bim = read.table("10k.bim")
library(data.table)
ld = fread("r2.ld")
ld = as.matrix(ld)
for(i in 1:20)
{
  idx = which(bim[,1] == i)
  m_ld = mean(ld[-idx, idx])
  print(c(i, m_ld))
}
```
# Two VCs
```sh
perl ../split_snps_to_2vcs.pl
perl ../make_grms.pl
```
# Simulate phenotypes
```{r }
hsq0 = 0.1
hsq1 = 0.4
nrep = 10
folder = "close"

library(data.table)

G = fread(paste(folder,"plink.0.rel",sep="/"))
G = as.matrix(G)
np = nrow(G)
y = t(chol(G)) * matrix(rnorm(np*nrep),np, nrep) * sqrt(hsq0)

G = fread(paste(folder,"plink.1.rel",sep="/"))
G = as.matrix(G)
y = y + t(chol(G)) * matrix(rnorm(np*nrep),np, nrep) * sqrt(hsq1)

y = y + matrix(rnorm(np*nrep),np, nrep) * sqrt(1-hsq0-hsq1)

id = read.table(paste(folder,"plink.0.rel.id",sep="/"))
out = matrix(nrow=np,ncol=(nrep+2))
out[,1:2] = id
out[,3:(2+nrep)] = y
colnames(out) = c("FID","IID",1:nrep)

write.table(out,file=paste(folder,"bolt.pheno.txt",sep="/"),quote=F,col.names=T,row.names=F,sep=" ")
write.table(out,file=paste(folder,"gcta.pheno.txt",sep="/"),quote=F,col.names=F,row.names=F,sep=" ")
write.table(out,file=paste(folder,"mph.pheno.csv",sep="/"),quote=F,col.names=F,row.names=F,sep=",")
```

