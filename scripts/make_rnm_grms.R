# Generate reaction-norm GRMs from a routine GRM and a quantitative covariate variable.
# Mar 7, 2024: initial release along with MPH v0.52.1

#!/usr/bin/env Rscript

source("mph_functs.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
    stop("Two arguments must be supplied:\n  input GRM filename prefix, input CSV of covariate, covariate name, highest polynomial order, and output filename prefix\n", call.=FALSE)
}

infile = args[1]
covarfile = args[2]
covarname = args[3]
order = as.numeric(args[4])
outfile = args[5]

covar = read.csv(covarfile, row.names = 1)
if(!(covarname %in% colnames(covar))) {
    stop(paste(covarname, "not found in the header of", covarfile ), call.=FALSE)
}

grm = read_grm(infile)
grm = grm + t(grm)
diag(grm) = diag(grm) / 2

common_iid <- intersect(rownames(grm), rownames(covar)[!is.na(covar[[covarname]])])
grm = grm[common_iid, common_iid]
x = covar[common_iid, covarname]

require(orthopolynom)
x = scaleX(x, u=-1, v=1)
# MTG2 does not scale this way. 
poly = legendre.polynomials(n=order, normalized=T)
p = matrix(unlist(polynomial.values(poly, x)), nrow=length(x))
for(i in 2:ncol(p)) {
	p[,i] = p[,i] / sqrt(mean(p[,i]^2))
	# p[,i] = x ^ (i-1)
	# MTG2 used the simple single-term polynomials.
}

for(i in 1:ncol(p)) {
    for(j in i:ncol(p)) {
        if(i == 1 && j == 1) {
            next
        } else if (i == j) {
           out_grm = ( p[,i] %*% t(p[,j]) ) * grm
        } else {
           out_grm = ( p[,i] %*% t(p[,j]) + p[,j] %*% t(p[,i]) ) * grm
        }
        write_grm(paste0(outfile, ".rnm.", paste(i-1,j-1,sep="x")), colnames(out_grm), out_grm)
        print(paste("GRM for order", i-1, "x", "order", j-1, "has been computed."))
    }
}

