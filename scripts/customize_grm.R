gbin = "test.grm.bin"
gid = "test.grm.id"

# Example data
# GRM
grm = matrix(1:25, ncol=5)
grm = grm + t(grm)
# IDs
id = c(1:5)
# sample size
np = 5
# number of markers set to 1 for pseudo GRMs
nm = 1

con = file(gbin, "wb")
writeBin(as.integer(np), con, size=4)
writeBin(as.numeric( c(nm,grm[lower.tri(grm, diag=TRUE)]) ), con, size=4)
close(con)

write.table(id, file=gid, quote=F, row.names=F, col.names=F)

# check .bin file
mph_file = file(gbin, "rb")
np = readBin(mph_file, n=1, what=integer(), size=4)
nm = readBin(mph_file, n=1, what=numeric(), size=4)
mph_G = matrix(0, nrow=np, ncol=np)
mph_G[lower.tri(mph_G, diag=TRUE)] = readBin(mph_file, n=np*(np+1)/2, what=numeric(), size=4)
close(mph_file)


