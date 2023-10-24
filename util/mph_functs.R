# Recalculate the estimates and SEs of PVEs and enrichments from MPH VC estimates
# vcfile: the .mq.vc.csv file produced by `mph --minque`.
# crossprod: the crossproduct of SNP incidence matrix and SNP weighting matrix. Its row names and columns should match the annotation categories of interest and the rows of vcfile, respectively.
# nsnps: the total number of SNPs. If `NA`, it is set to the first GRM's number of SNPs in vcfile.
# annot.size: a list of the number of SNPs in each annotation category listed in the row names of crossprod. It can be computed from the column-wise sum of the corresponding SNP incidence matrix. If `NA`, it is set to the `m` column of vcfile.
recalculate_enrichments <- function(vcfile, crossprod, nsnps=NA, annot.size=NA) {
    if(is.na(vcfile)) {
        stop("vcfile must be specifiled.")
    }
    if(length(crossprod) == 1 && is.na(crossprod)) {
        stop("crossprod must be specified.")
    }

    mq = read.csv(vcfile)
    mq = mq[-nrow(mq),]
    if(ncol(crossprod) > nrow(mq)) {
        stop(paste("ncol(crossprod) should equal the number of genetic VCs in", vcfile))
    }
    mq = mq[1:ncol(crossprod),]

    # PVE estimate
    var = crossprod %*% as.matrix(mq[,9:(8+nrow(mq))]) %*% t(crossprod)
    pve.est = crossprod %*% mq$enrichment / sum(mq$m)
    pve.se = sqrt(diag(var)) / sum(mq$m) 

    # enrichment estimate
    if(is.na(nsnps)) {
        nsnps = mq$m[1]
    }
    if(length(annot.size) == 1 && is.na(annot.size)) {
        annot.size = mq$m
    }
    if(length(annot.size) != nrow(crossprod)) {
        stop("annot.size must correspond to the rows of crossprod.")
    }
    e.est = pve.est / (annot.size/nsnps)
    e.se = pve.se / (annot.size/nsnps)

    out = cbind(annot.size/nsnps, pve.est, pve.se, e.est, e.se)
    rownames(out) = rownames(crossprod)
    colnames(out) = c("prop", "pve", "pve.se", "enrichment", "enrichment.se")
    return(out)
}

# Write GRM into MPH format files (.iid and .bin)
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

# Read MPH GRM files (.iid and .bin)
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
