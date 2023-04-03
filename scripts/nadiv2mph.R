library(nadiv)

# individual sire dam
pedfile = "../renadd01.ped"

ped = read.table(pedfile)
ped = ped[,c(1,3,2)]

goodped = prepPed(ped)
mat = makeD(goodped, parallel = TRUE, ncores = getOption("mc.cores", 10L), invertD = FALSE, returnA = FALSE, det = FALSE)

grm = as.matrix(mat$D)
iid = colnames(grm)
write_grm("dom.d", iid, grm)
