#!/usr/bin/env Rscript
source("grm_io.R")
library(data.table)

args = commandArgs(trailingOnly=TRUE)
# test if there are two arguments: if not, return an error
if (length(args) != 2) {
  stop("Two arguments must be supplied:  \ninput txt filename and output filename prefix\n")
}

txtfile = args[1]
prefix = args[2]
print(paste("Output filename prefix:", prefix))

grm = as.matrix(fread(txtfile, head=F, check.names=F), rownames=1)
print(paste("Completed reading the input file:", txtfile))

iid = rownames(grm)

write_grm(prefix, iid, grm)
