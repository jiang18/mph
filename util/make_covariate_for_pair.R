# Generates a new covariate file for a trait pair
# Oct 27, 2023: initial release along with MPH v0.49.2

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    stop("Two arguments must be supplied:\n  input covariate filename and output filename\n", call.=FALSE)
}

infile = args[1]
outfile = args[2]

covar = read.csv(infile)
nr = nrow(covar)
nc = ncol(covar)-1

iid = c(paste0("t1.", covar[,1]), paste0("t2.", covar[,1]))
out = cbind(iid, matrix(0, nrow=2*nr, ncol=2*nc))
out[1:nr, 2:(nc+1)] = as.matrix(covar[,-1])
out[(nr+1):(2*nr), (nc+2):(2*nc+1)] = as.matrix(covar[,-1])

colnames(out) = c("IID", paste0("t1.", colnames(covar)[-1]), paste0("t2.", colnames(covar)[-1]))
write.csv(out, file=outfile, quote=F, row.names=F)
