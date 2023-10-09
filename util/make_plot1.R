d = read.csv("../11-100.parse.txt")

tt = d$trt[which(d$diff<0.01)]
tt = c(1:10, tt)

summary = matrix(nrow=81, ncol=25)
j = 1
for (i in tt) {
	dd = read.csv( paste("../50k/hsq90",i,"200.1.mq.csv", sep=".") )
	summary[j,] = dd$var[-26]
	j = j+1
}


sim_true = read.table("../ng2015/sim_true.txt",sep="\t",quote="")
dt = summary
colnames(dt) = unlist(strsplit(sim_true[,1], "[.]"))[c(1,seq(2,48,by=2))]

par(mar=c(2,12,1,1))
boxplot(dt, horizontal=TRUE, las=1)
points(colMeans(dt), 1:ncol(dt), col = "red", pch=1)
points(sim_true[,2], 1:ncol(dt), col = "blue", pch=3)
abline(v=1,col="blue")
legend("topright", legend=c("Estimate mean", "True value"), col=c("red", "blue"), pch=c(1,3))

dev.off()



###########################


library(data.table)
sw = fread("../funct.snp_info.csv", head=T, sep=",")

rname = sw$snp
sw = sw[,-1]
sw[is.na(sw)] = 0
sw = as.matrix(sw)

main = which( !(grepl("500", colnames(sw)) | grepl("peak", colnames(sw))) )
sw = sw[,main]

rownames(sw) = rname

# MPH output
mq = read.csv("../50k/hsq90.1.200.1.mq.csv")
mq = mq[-nrow(mq),]

# heritability estimate
lhs = t(sw) %*% sw
var = lhs %*% as.matrix(mq[,9:33]) %*% t(lhs)
h.est = lhs %*% mq$enrichment / sum(mq$m)
h.se = sqrt(diag(var)) / sum(mq$m) 

# enrichment estimate
e.est = h.est / (mq$m/mq$m[1])
e.se = h.se / (mq$m/mq$m[1])

ee = cbind(e.est[-1], e.se[-1])
rownames(ee) = unlist(strsplit(rownames(ee), "[.]"))[seq(1,48,by=2)]

df = data.frame(FunctionalAnnotation=rownames(ee), Enrichment=ee[,1], SE=ee[,2])


p<- ggplot(df, aes(x=FunctionalAnnotation, y=Enrichment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Enrichment-SE, ymax=Enrichment+SE), width=.2,
                 position=position_dodge(.9)) 
p<- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p<- p+ geom_hline(yintercept=1, linetype="dashed", color = "red") + xlab("")
print(p)
dev.off()
