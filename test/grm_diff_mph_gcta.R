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
  closeAllConnections()
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}
##########################################################

bin = ReadGRMBin(paste(folder,"gcta.0",sep="/"))
np = length(bin$diag)
G = matrix(0, nrow=np, ncol=np)
G[upper.tri(G)] = bin$off
diag(G) = bin$diag

mph_file = file(paste(folder,"mph.0",sep="/"), "rb")
readBin(mph_file, n=1, what=integer(), size=4)
readBin(mph_file, n=1, what=numeric(), size=4)
mph_G = matrix(0, nrow=np, ncol=np)
mph_G[lower.tri(G, diag=TRUE)] = readBin(mph_file, n=np*(np+1)/2, what=numeric(), size=4)
mph_G = t(mph_G)

max(abs(mph_G - G))

