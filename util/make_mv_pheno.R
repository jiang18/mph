# Generate multivariate phenotype and covariate files
# Oct 31, 2023: initial release along with MPH v0.49.2

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("At least four arguments must be supplied:\n  input phenotype CSV filename, 2+ trait names, and output CSV filename prefix\n", call.=FALSE)
}

infile = args[1]
mvtrt = args[2:(length(args)-1)]
outfile = args[length(args)]

if(any(duplicated(mvtrt))) {
    stop("Trait names must differ from each other.\n", call.=FALSE)
}

phe = read.csv(infile, row.names=1, check.names = FALSE)
trt = colnames(phe)
if(! all(mvtrt %in% trt) ) {
    stop(paste("One or more trait names are unavailable in", infile), call.=FALSE)
}

phe = na.omit( phe[, mvtrt] )
if(nrow(phe) == 0) {
    stop(paste("No observations are available for all traits:", mvtrt), call.=FALSE)
}

out = cbind(apply(phe, 2, mean), apply(phe, 2, sd))
colnames(out) = c("mean", "std")
write.csv(out, file=paste(outfile, paste(mvtrt, collapse = "."), "std.csv", sep="."), quote=F, row.names=T)

dim = nrow(phe)
ntrt = length(mvtrt)
iid = paste(rep(1:ntrt, each=dim), rep(rownames(phe), ntrt), sep=".")

out = cbind(iid, as.vector(apply(phe, 2, scale)))
colnames(out) = c("IID", "scaled")
write.csv(out, file=paste(outfile, paste(mvtrt, collapse = "."), "pheno.csv", sep="."), quote=F, row.names=F)

out = matrix(0, nrow=dim*ntrt, ncol=ntrt+1)
out[, 1] = iid
for(i in 1:ntrt) {
    out[(1+(i-1)*dim):(i*dim), i+1] = 1
}
colnames(out) = c("IID", mvtrt)
write.csv(out, file=paste(outfile, paste(mvtrt, collapse = "."), "covar.csv", sep="."), quote=F, row.names=F, na="")
