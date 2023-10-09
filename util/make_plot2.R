d = read.csv("Milk.mq.csv")
d = d[-nrow(d),]

df = data.frame(MAFbin=c("0.005-0.05","0.05-0.10","0.10-0.15","0.15-0.20","0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40","0.40-0.45","0.45-0.50"), Enrichment=d$enrichment, SE=d$seE)


p<- ggplot(df, aes(x=MAFbin, y=Enrichment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Enrichment-SE, ymax=Enrichment+SE), width=.2,
                 position=position_dodge(.9)) 
p<- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p<- p+ geom_hline(yintercept=1, linetype="dashed", color = "red") + xlab("MAF bin") + ylab("Heritability enrichment")
print(p)
dev.off()

###################
require(gridExtra)

d = read.csv("Fat_Percent.mq.csv")
d = d[-nrow(d),]

cum = rep(0,5)
cum_se = rep(0,5)
for(i in seq(2,10,by=2)) {
	cum[i/2] = sum(d$m[1:i] * d$enrichment[1:i])/sum(d$m)
	temp = d$m[1:i] / sum(d$m)
	cum_se[i/2] = sqrt(sum(temp * ( as.matrix(d[1:i,9:(8+i)]) %*% temp )) )
}
df = data.frame(MAF=c(1:5)/10, Cumulative=cum, SE=cum_se)
p<- ggplot(df, aes(x=MAF, y=Cumulative)) + 
  geom_point() + xlim(0,0.505) + ylim(0,1.01) + 
  geom_errorbar(aes(ymin=Cumulative, ymax=Cumulative+SE), width=.01,
                 position=position_dodge(.9)) 
p<- p+ geom_segment(aes(x=0,xend=0.5,y=0,yend=1), color = "red") + xlab("MAF") + ylab("Cumulative contribution to\n genetic variance")
plot1 = p


cum = rep(0,5)
cum_se = rep(0,5)
for(i in seq(2,10,by=2)) {
	cum[i/2] = sum(d$m[(i-1):i] * d$enrichment[(i-1):i])/sum(d$m)
	temp = d$m[(i-1):i] / sum(d$m)
	cum_se[i/2] = sqrt(sum(temp * ( as.matrix(d[(i-1):i,(7+i):(8+i)]) %*% temp )) )
}
df = data.frame(MAFbin=c("0.005-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5"), PVE=cum, SE=cum_se)
p<- ggplot(df, aes(x=MAFbin, y=PVE)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=PVE, ymax=PVE+SE), width=.2,
                 position=position_dodge(.9)) 
p<- p+ xlab("MAF bin") + ylab("Contribution to genetic variance")
plot2 = p

grid.arrange(plot1, plot2, ncol=2, heights=unit(4, "in"))
dev.off()


#####


d = read.csv("Pro_Percent.mq.csv")
d = d[-nrow(d),]

cum = rep(0,10)
cum_se = rep(0,10)
for(i in seq(1,10,by=1)) {
	cum[i] = sum(d$m[1:i] * d$enrichment[1:i])/sum(d$m)
	temp = d$m[1:i] / sum(d$m)
	cum_se[i] = sqrt(sum(temp * ( as.matrix(d[1:i,9:(8+i)]) %*% temp )) )
}
df = data.frame(MAF=c(1:10)/20, Cumulative=cum, SE=cum_se)
p<- ggplot(df, aes(x=MAF, y=Cumulative)) + 
  geom_line() + geom_point() + xlim(0,0.505) + ylim(0,1.01) + 
  geom_errorbar(aes(ymin=Cumulative-SE, ymax=Cumulative+SE), width=.01,
                 position=position_dodge(.9)) 
p<- p+ geom_segment(aes(x=0,xend=0.5,y=0,yend=1), color = "red") + xlab("MAF") + ylab("Cumulative contribution to\n genetic variance")
plot1 = p


df = data.frame(MAFbin=c("0.005-0.05","0.05-0.10","0.10-0.15","0.15-0.20","0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40","0.40-0.45","0.45-0.50"), PVE=d$m * d$enrichment / sum(d$m), SE=d$seE * d$m / sum(d$m))
p<- ggplot(df, aes(x=MAFbin, y=PVE)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=PVE, ymax=PVE+SE), width=.2,
                 position=position_dodge(.9)) 
p<- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p<- p+ xlab("MAF bin") + ylab("Contribution to genetic variance")
plot2 = p

grid.arrange(plot1, plot2, ncol=2, heights=unit(4, "in"))
dev.off()



