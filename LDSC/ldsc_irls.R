args <- commandArgs(trailingOnly=TRUE)
in_prefix = args[1]

ldsc_file = "close/ldsc"
assoc_file = paste(in_prefix,".qassoc",sep="")
out_file = paste(in_prefix,".irls",sep="")

pow = -1
n = 10000
ldsc = read.table(ldsc_file, head=T)
assoc = read.table(assoc_file, head=T)
ldsc = as.matrix(ldsc)

idx_na = which(is.na(assoc$T))
assoc = assoc[-idx_na,]

output = matrix(nrow=3, ncol=ncol(ldsc)+1)
colnames(output) = c(colnames(ldsc), "C")
rownames(output) = c("m","unitC","varC")

m = ldsc[1,]
ldsc = ldsc[-1,]
ldsc = ldsc[-idx_na,]

output[1,] = c(m,sum(m))

# genome-wide ld score
gwsc = rowSums(ldsc)
chisq = assoc$T **2

pw = gwsc **pow
pw = pw/sum(pw)

ll = t(t(ldsc)*n-m)
irls_unit <- function(theta) {
	theta_old = theta
	for(i in 1:10) {
		ltheta = ll %*% theta_old + 1
		weight = pw*(2*chisq/ltheta - 1)/ltheta/ltheta
		delta = pw*(chisq - ltheta)/ltheta/ltheta
		z = (ltheta - 1) + delta/weight
		xvz = t(ll) %*% (weight[,1] * z)
		xvx = t(ll) %*% (weight[,1] * ll)
		
		theta_new = solve(xvx, xvz)
		theta_old = theta_new
		# print(theta_new)
	}
	return(theta_new)
}
theta = c(1e-3/m)
theta = irls_unit(theta)
output[2,] = c(theta,1)

ll = cbind(ll,1)
irls <- function(theta) {
	theta_old = theta
	for(i in 1:10) {
		ltheta = ll %*% theta_old
		weight = pw*(2*chisq/ltheta - 1)/ltheta/ltheta
		delta = pw*(chisq - ltheta)/ltheta/ltheta
		z = ltheta + delta/weight
		xvz = t(ll) %*% (weight[,1] * z)
		xvx = t(ll) %*% (weight[,1] * ll)
		
		theta_new = solve(xvx, xvz)
		theta_old = theta_new
		# print(theta_new)
	}
	return(theta_new)
}

theta = c(1e-3/m,1)
theta = irls(theta)
output[3,] = theta

output[2,1:2] = output[1,1:2] * output[2,1:2]
output[3,1:2] = output[1,1:2] * output[3,1:2]

write.table(output, file=out_file, quote=F)



