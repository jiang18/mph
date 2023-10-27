# Generate new phenotype and covariate files for a trait pair
# Oct 27, 2023: initial release along with MPH v0.49.2

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    stop("Four arguments must be supplied:\n  input phenotype filename, trait 1, trait 2, and output filename prefix\n", call.=FALSE)
}

infile = args[1]
t1 = args[2]
t2 = args[3]
outfile = args[4]

if(t1 == t2) {
    stop("The two trait names are identical.\n", call.=FALSE)
}

phe = read.csv(infile)
trt = colnames(phe)
if(!(t1 %in% trt)) {
    stop(paste(t1, "is not available in", infile), call.=FALSE)
}
if(!(t2 %in% trt)) {
    stop(paste(t2, "is not available in", infile), call.=FALSE)
}

phe = phe[!is.na(phe[[t1]]) & !is.na(phe[[t2]]), ]
if(nrow(phe) == 0) {
    stop(paste("No observations are available for", t1, t2), call.=FALSE)
}

out = cbind(c(t1,t2), c(sd(phe[[t1]]), sd(phe[[t2]])))
colnames(out) = c("trait", "std")
write.csv(out, file=paste(outfile, t1, t2, "std.csv", sep="."), quote=F, row.names=F)

iid = c(paste0("t1.", phe[,1]), paste0("t2.", phe[,1]))
out = cbind(iid, c(scale(phe[[t1]]), scale(phe[[t2]])))
colnames(out) = c("IID", "scaled")
write.csv(out, file=paste(outfile, t1, t2, "pheno.csv", sep="."), quote=F, row.names=F)

out = cbind(iid, c(rep(1, nrow(phe)), rep(0, nrow(phe))))
out = cbind(out, c(rep(0, nrow(phe)), rep(1, nrow(phe))))
colnames(out) = c("IID", "t1", "t2")
write.csv(out, file=paste(outfile, t1, t2, "covar.csv", sep="."), quote=F, row.names=F)
