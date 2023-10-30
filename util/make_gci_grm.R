# Generate a genotype–covariate interaction GRM from a routine GRM and an incidence matrix
# Oct 29, 2023: initial release along with MPH v0.49.2

#!/usr/bin/env Rscript

source("mph_functs.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    stop("Two arguments must be supplied:\n  input GRM filename prefix, input CSV of covariate incidence, and output filename prefix\n", call.=FALSE)
}

infile = args[1]
covarfile = args[2]
outfile = args[3]

covar = read.csv(covarfile, row.names = 1)
if( ! all(rowSums(covar) == 1) ) {
    stop(paste(covarfile, "does not seem to be a CSV file of covariate incidence.\n"))
}

grm = read_grm(infile)
grm = grm + t(grm)
diag(grm) = diag(grm) / 2

common_iid <- intersect(rownames(grm), rownames(covar))
grm = grm[common_iid, common_iid]
covar = as.matrix( covar[common_iid, ] )
grm = grm * tcrossprod(covar)

iid = rownames(grm)
write_grm(paste0(outfile, ".gci"), iid, grm)
