# R script to write GRM into MPH files (.iid and .bin)
write_grm=function(prefix, iid, grm){
  idfile = paste0(prefix, ".grm.iid")
  binfile = paste0(prefix, ".grm.bin")
  if(length(iid) != ncol(grm)) {
    print("Error: the size of iid is not equal to ncol of grm.")
    return(1)
  }
  write.table(iid, file=idfile, quote=F, row.names=F, col.names=F)
  print(paste0("Completed writing the iid file: ", prefix, ".grm.iid"))
  
  # sample size
  np = length(iid)
  # number of markers set to 1 for pseudo GRMs
  nm = 1
  con = file(binfile, "wb")
  writeBin(as.integer(np), con, size=4)
  writeBin(as.numeric( c(nm,grm[lower.tri(grm, diag=TRUE)]) ), con, size=4)
  close(con)
  print(paste0("Completed writing the bin file: ", prefix, ".grm.bin"))
}

# R script to read MPH GRM files (.iid and .bin)
read_grm=function(prefix){
  idfile = paste0(prefix, ".grm.iid")
  binfile = paste0(prefix, ".grm.bin")
  iid = read.table(idfile)
  iid = iid[,1]
  mph_file = file(binfile, "rb")
  np = readBin(mph_file, n=1, what=integer(), size=4)
  if(length(iid) != np) {
    print("Error: the size of iid is not equal to ncol of grm.")
    return(1)
  }
  nm = readBin(mph_file, n=1, what=numeric(), size=4)
  grm = matrix(0, nrow=np, ncol=np)
  grm[lower.tri(grm, diag=TRUE)] = readBin(mph_file, n=np*(np+1)/2, what=numeric(), size=4)
  close(mph_file)
  
  grm = grm / nm
  colnames(grm) = iid
  rownames(grm) = iid

  return(grm)
}
