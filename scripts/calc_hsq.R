library(data.table)
sw = fread("funct.snp_info.csv", head=T, sep=",")

rname = sw$snp
sw = sw[,-1]
sw[is.na(sw)] = 0
sw = as.matrix(sw)

main = which( !(grepl("500", colnames(sw)) | grepl("peak", colnames(sw))) )
sw = sw[,main]

rownames(sw) = rname

ng = read.table("ng2015/enrichment.txt",sep="\t",head=T)
rhs = c(1,ng[,3])
lhs = t(sw) %*% sw
fold = solve(lhs, rhs)
fold = colSums(sw) * fold

fold[1] = fold[1] + 0.3
snpvar = sw %*% (fold/colSums(sw))
summary(snpvar)
prop = t(sw) %*% snpvar / sum(fold)
cbind(prop, c(1,ng[,3]))

write.table(fold, file="ng2015/sim_true.txt", col.names=F, row.names=T, quote=F, sep="\t")


####################################

# true values
ng = read.table("../ng2015/enrichment.txt",sep="\t",head=T)
true = c(1,ng[,3])

# W matrix
sw = fread("../funct.snp_info.csv", head=T, sep=",")
rname = sw$snp
sw = sw[,-1]
sw[is.na(sw)] = 0
sw = as.matrix(sw)
main = which( !(grepl("500", colnames(sw)) | grepl("peak", colnames(sw))) )
sw = sw[,main]
rownames(sw) = rname

# MPH output
mq = read.csv("hsq90.1.200.mq.csv")
mq = mq[-nrow(mq),]

# heritability estimate
lhs = t(sw) %*% sw
var = lhs %*% as.matrix(mq[,9:33]) %*% t(lhs)
h.est = lhs %*% mq$enrichment / sum(mq$m)
h.se = sqrt(diag(var)) / sum(mq$m) 

# enrichment estimate
e.est = h.est / (mq$m/mq$m[1])
e.se = h.se / (mq$m/mq$m[1])

cbind(ng, h.est[-1], h.se[-1], e.est[-1], e.se[-1])

#######################################

par(mfrow=c(2,1))
plot(ng[,3], h.est[-1])
abline(a=0, b=1)
plot(ng[,4], h.se[-1])
abline(a=0, b=1)
dev.off()


