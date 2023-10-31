# Generate multivariate genomic relationship matrices
# Oct 31, 2023: initial release along with MPH v0.49.2

#!/usr/bin/env Rscript

source("mph_functs.R")

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
    stop("Two or three arguments must be supplied:\n  input GRM filename prefix, output filename prefix, and number of traits (optional).\n", call.=FALSE)
}

if(ntrt > 9) {
    stop("Too many traits: The max number of traits is limited to 9.\n")
}

grm = read_grm(infile)
grm = grm + t(grm)
diag(grm) = diag(grm) / 2

dim = dim(grm)[1]
iid = paste(rep(1:ntrt, each=dim), rep(colnames(grm), ntrt), sep=".")

for(i in 1:ntrt) {
    full = matrix(0, nrow=dim*ntrt, ncol=dim*ntrt)
    full[(1+(i-1)*dim):(dim*i), (1+(i-1)*dim):(dim*i)] = grm
    write_grm(paste0(outfile, ".G", i, i), iid, full)
}

for(i in 1:(ntrt-1)) {
    for(j in (i+1):ntrt) {
        full = matrix(0, nrow=dim*ntrt, ncol=dim*ntrt)
        full[(1+(i-1)*dim):(dim*i), (1+(j-1)*dim):(dim*j)] = grm
        full[(1+(j-1)*dim):(dim*j), (1+(i-1)*dim):(dim*i)] = grm
        write_grm(paste0(outfile, ".G", i, j), iid, full)
    }
}
