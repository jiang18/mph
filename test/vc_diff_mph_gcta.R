mph_out = "10.k.30.mq.csv"
gcta_out = "10.k.30.hsq"
nchr = 30

y = read.csv( mph_out )
x = read.table(gcta_out,head=T,nrow=nchr+1)
y = y[-nrow(y),]
x = x[-nrow(x),]

if(FALSE) {
  var = read.table( paste(gcta_out,".var",sep="") )
  var = as.matrix(var)
  var = var[-nrow(var),-ncol(var)]
  m = y$m
  h = x$Variance
  one = rep(1,nrow(var))
  rho = sum(m)/sum(h)*(h/m)
  J = sum(m)/sum(h)**2 * ( sum(h)*diag(1/m) - diag(1/m) %*% h %*% t(one) )
  se_rho = sqrt(diag(J %*% var %*% t(J)))
}

par(mfrow=c(3,1))
summary(lm(y$var~x$Variance))
cor(y$var, x$Variance)
plot(x$Variance, y$var)
abline(a=0,b=1,col="blue")

summary(lm(y$seV~x$SE))
cor(y$seV, x$SE)
plot(x$SE, y$seV)
abline(a=0,b=1,col="blue")

tstat.y = y$var/y$seV
tstat.x = x$Variance/x$SE
summary(lm(tstat.y~tstat.x))
cor(tstat.y, tstat.x)
plot(tstat.x, tstat.y)
abline(a=0,b=1,col="blue")
dev.off()
