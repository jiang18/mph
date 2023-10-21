#!/usr/bin/env Rscript
source("grm_io.R")
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
    stop("Two arguments must be supplied:\n  input GRM filename prefix and output filename prefix\n", call.=FALSE)
}

infile = args[1]
outfile = args[2]

keep = read.table(paste0(infile, ".grm.iid"))
dim = length(keep[,1])
iid = c(paste0("t1.", keep[,1]), paste0("t2.", keep[,1]))

full = matrix(0, nrow=dim*2, ncol=dim*2)
diag(full[1:dim, 1:dim]) = 1
write_grm(paste0(outfile, ".E1d"), iid, full)

full = matrix(0, nrow=dim*2, ncol=dim*2)
diag(full[1:dim, (1+dim):(dim*2)]) = 1
diag(full[(1+dim):(dim*2), 1:dim]) = 1
write_grm(paste0(outfile, ".E12"), iid, full)
