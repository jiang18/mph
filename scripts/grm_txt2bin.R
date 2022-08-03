#!/usr/bin/env Rscript

library(data.table)

args = commandArgs(trailingOnly=TRUE)
# test if there are two arguments: if not, return an error
if (length(args) != 2) {
  stop("Two arguments must be supplied (input and output).n", call.=FALSE)
}

txtfile = args[1]
prefix = args[2]
print("Input filename prefix:", prefix)

idfile = paste0(prefix, ".indi")
binfile = paste0(prefix, ".bin")

grm = as.matrix(fread(txtfile, head=F, check.names=F), rownames=1)
id = rownames(grm)
write.table(id, file=idfile, quote=F, row.names=F, col.names=F)

# sample size
np = length(id)
# number of markers set to 1 for pseudo GRMs
nm = 1
con = file(binfile, "wb")
writeBin(as.integer(np), con, size=4)
writeBin(as.numeric( c(nm,grm[lower.tri(grm, diag=TRUE)]) ), con, size=4)
close(con)
