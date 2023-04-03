# use AGHmatrix to construct numerator relationship matrix with renumf90 ped file

library(AGHmatrix)

args = commandArgs(trailingOnly=TRUE)
# test if there are three arguments: if not, return an error
if (length(args) != 3) {
  stop("Three arguments must be supplied: renumf90 ped filename, mph grm prefix, and A (or D).\n")
}

pedfile = args[1]
prefix = args[2]
pedmat = args[3]
if (pedmat != "A" && pedmat != "D") {
  stop("The 3rd argument must be A or D.\n")
}
print(paste("Output filename prefix:", prefix))

ped = read.table(pedfile)
print("Completed reading the input file.")
ped = ped[,1:3]
if (pedmat == "A") {
  grm = Amatrix(ped, ploidy=2)
} else {
  grm = Amatrix(ped, ploidy=2, dominance=TRUE)
}

idfile = paste0(prefix, ".grm.iid")
binfile = paste0(prefix, ".grm.bin")

id = colnames(grm)
write.table(id, file=idfile, quote=F, row.names=F, col.names=F)
print("Completed writing the iid file.")

# sample size
np = length(id)
# number of markers set to 1 for pseudo GRMs
nm = 1
con = file(binfile, "wb")
writeBin(as.integer(np), con, size=4)
writeBin(as.numeric( c(nm,grm[lower.tri(grm, diag=TRUE)]) ), con, size=4)
close(con)
print("Completed writing the bin file.")
