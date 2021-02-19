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
cp ../chips.txt ./
cp ../markersim.options ./
cp ../genosim.options. ./
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
# Copied from https://cnsgenomics.com/software/gcta/#MakingaGRM
# R script to read the GRM binary file
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}
##########################################################

hsq0 = 0.1
hsq1 = 0.4
nrep = 10
folder = "close"

bin = ReadGRMBin(paste(folder,"gcta.0",sep="/"))
np = length(bin$diag)
G = matrix(0, nrow=np, ncol=np)
G[upper.tri(G)] = bin$off
G = G + t(G)
diag(G) = bin$diag + 1e-6
y = t(chol(G)) %*% matrix(rnorm(np*nrep),np, nrep) * sqrt(hsq0)

bin = ReadGRMBin(paste(folder,"gcta.1",sep="/"))
G[,] = 0
G[upper.tri(G)] = bin$off
G = G + t(G)
diag(G) = bin$diag + 1e-6
y = y + t(chol(G)) %*% matrix(rnorm(np*nrep),np, nrep) * sqrt(hsq1)

y = y + matrix(rnorm(np*nrep),np, nrep) * sqrt(1-hsq0-hsq1)

closeAllConnections()

out = matrix(nrow=np,ncol=(nrep+2))
out[,1:2] = as.matrix(bin$id)
out[,3:(2+nrep)] = y
colnames(out) = c("FID","IID",1:nrep)

write.table(out,file=paste(folder,"bolt.pheno.txt",sep="/"),quote=F,col.names=T,row.names=F,sep=" ")
write.table(out,file=paste(folder,"gcta.pheno.txt",sep="/"),quote=F,col.names=F,row.names=F,sep=" ")
out = out[,-1]
write.table(out,file=paste(folder,"mph.pheno.csv",sep="/"),quote=F,col.names=T,row.names=F,sep=",")
```
# Estimate VCs
```sh
perl ../estimate_vcs.pl
Rscript --no-save ../cal_scores.R
perl ../run_ldsc.pl
```
