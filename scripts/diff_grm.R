# diff grms

grm1 = "../genome";
grm2 = "Pre_LNC.txt"
out = "pre_remain"

id1 = read.table(paste(grm1,"grm.indi", sep="."))
id2 = read.table(paste(grm2,"grm.indi", sep="."))

sum(id1[,1] == id2[,1])


gbin = paste(grm1,"grm.bin", sep=".")
mph_file = file(gbin, "rb")
np = readBin(mph_file, n=1, what=integer(), size=4)
nm1 = readBin(mph_file, n=1, what=numeric(), size=4)
mph_G1 = matrix(0, nrow=np, ncol=np)
mph_G1[lower.tri(mph_G1, diag=TRUE)] = readBin(mph_file, n=np*(np+1)/2, what=numeric(), size=4)
close(mph_file)

gbin = paste(grm2,"grm.bin", sep=".")
mph_file = file(gbin, "rb")
np = readBin(mph_file, n=1, what=integer(), size=4)
nm2 = readBin(mph_file, n=1, what=numeric(), size=4)
mph_G2 = matrix(0, nrow=np, ncol=np)
mph_G2[lower.tri(mph_G2, diag=TRUE)] = readBin(mph_file, n=np*(np+1)/2, what=numeric(), size=4)
close(mph_file)

grm = mph_G1 - mph_G2
nm = nm1 - nm2
gbin = paste(out,"grm.bin", sep=".")
con = file(gbin, "wb")
writeBin(as.integer(np), con, size=4)
writeBin(as.numeric( c(nm,grm[lower.tri(grm, diag=TRUE)]) ), con, size=4)
close(con)
