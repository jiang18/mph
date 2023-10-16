#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
    stop("Two arguments must be supplied:\n  input GRM filename prefix and output filename prefix\n", call.=FALSE)
}

infile = args[1]
outfile = args[2]

source("grm_io.R")

grm = read_grm(infile)
grm = grm + t(grm)
diag(grm) = diag(grm) / 2

dim = dim(grm)[1]
iid = c(paste0("t1.", colnames(grm)), paste0("t2.", colnames(grm)))
full = matrix(0, nrow=dim*2, ncol=dim*2)
full[1:dim, 1:dim] = grm
write_grm(paste0(outfile, ".G.t1"), iid, full)

full = matrix(0, nrow=dim*2, ncol=dim*2)
full[(1+dim):(dim*2), (1+dim):(dim*2)] = grm
write_grm(paste0(outfile, ".G.t2"), iid, full)

full = matrix(0, nrow=dim*2, ncol=dim*2)
full[1:dim, (1+dim):(dim*2)] = grm
full[(1+dim):(dim*2), 1:dim] = grm
write_grm(paste0(outfile, ".G.t1t2"), iid, full)

full = matrix(0, nrow=dim*2, ncol=dim*2)
diag(full[1:dim, 1:dim]) = 1
write_grm(paste0(outfile, ".E.t1"), iid, full)

full = matrix(0, nrow=dim*2, ncol=dim*2)
diag(full[1:dim, (1+dim):(dim*2)]) = 1
diag(full[(1+dim):(dim*2), 1:dim]) = 1
write_grm(paste0(outfile, ".E.t1t2"), iid, full)
