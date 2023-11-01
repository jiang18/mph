# Generates a multivariate covariate file
# Oct 31, 2023: initial release along with MPH v0.49.2

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 2) {
    infile = args[1]
    outfile = args[2]
    ntrt = 2
} else if(length(args) == 3) {
    infile = args[1]
    outfile = args[2]
    ntrt = as.integer(args[3])
} else {
    stop("Two or three arguments must be supplied:\n  input covariate CSV filename, output CSV filename, and number of traits (optional)\n", call.=FALSE)
}

covar = read.csv(infile, check.names = FALSE)
dim = nrow(covar)
nc = ncol(covar)-1
iid = paste(rep(1:ntrt, each=dim), rep(covar[,1], ntrt), sep=".")

out = matrix(0, nrow=dim*ntrt, ncol=ntrt*nc+1)
out[,1] = iid
for(i in 1:ntrt) {
    out[(1+(i-1)*dim):(i*dim), (2+(i-1)*nc):(1+i*nc)] = as.matrix(covar[,-1])
}
colnames(out) = c("IID", paste(rep(1:ntrt, each=nc), rep(colnames(covar)[-1], ntrt), sep="."))
write.csv(out, file=outfile, quote=F, row.names=F, na="")
