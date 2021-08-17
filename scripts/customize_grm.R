gbin = "test.grm.bin"
gid = "test.grm.id"

# Example data
# GRM
grm = matrix(1:25, ncol=5)
grm = grm + t(grm)
# IDs
id = c(1:5)
# sample size
np = as.integer(5)
# number of markers set to 1 for pseudo GRMs
nm = 1

con = file(gbin, "wb")
writeBin(con=gbin, object=np)
writeBin(con=gbin, object=as.numeric( c(nm,grm[lower.tri(grm, diag=TRUE)]) ), size=4)
closeAllConnections()

write.table(id, file=gid, quote=F, row.names=F, col.names=F)

# check files

