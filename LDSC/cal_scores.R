folder = "close"

library(data.table)
ld = fread("r2.ld")
ld = as.matrix(ld)

snp = read.table(paste(folder,"bolt.modelSnps.txt", sep="/"))
weight = unique(snp[,2])
scores = matrix(nrow=nrow(ld),ncol=length(weight))
colnames(scores) = weight
num = rep(0,length(weight))
for(j in 1:length(weight)) {
	idx = which(snp[,2] == weight[j])
        num_na = sum(is.na(ld[1,idx]))
        num[j] = length(idx) - num_na
}
ld[is.na(ld)] = 0
for(j in 1:length(weight)) {
	idx = which(snp[,2] == weight[j])
	scores[,j] = rowSums(ld[,idx])
}
scores = rbind(num,scores)
write.table(scores, file=paste(folder,"ldsc",sep="/"), quote=F, row.names=F, col.names=T)




