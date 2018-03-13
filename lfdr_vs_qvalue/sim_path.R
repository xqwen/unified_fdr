################## locfdr vs qvalue ##################
source("ufdr.R")
source("mylocfdr.R")


##### simulated data, ZA ###

N = 20000
pi0 = 0.55
n1 = round((1-pi0)*N)
n2 = N-n1
tv = c(rep(1,n1),rep(0,n2))


shape = 1.0 
scale = 10
set.seed(123)

z1 = sqrt(rgamma(n1, shape=shape, scale = scale))*sample(c(-1,1),n1,replace=T)
z2 = rnorm(n2)
zval = c(z1,z2)
pval = 2*pnorm(-abs(zval))


hpi0 = qvalue(pval)$pi0
lfdr.rst = mylocfdr(zval,pi0=hpi0)

path_lfdr = get_rej_path_lfdr(lfdr.rst$fdr)
path_qval = get_rej_path_pval(pval,hpi0)

pdf(file="lfdr_vs_qval_k_1.pdf", width=5, height=5.5,bg="white")
plot(path_lfdr ~ path_qval, pch=16,cex=0.5, xlab="rejection path of q-value", ylab = "rejection path of local fdr", main = "k = 1.0")
abline(0,1,lty=2, col="darkgreen")
dev.off()






# heavy tailed distribution, UA

N = 20000
pi0 = 0.55
n1 = round((1-pi0)*N)
n2 = N-n1
tv = c(rep(1,n1),rep(0,n2))


shape = 0.3
scale = 10
set.seed(123)

z1 = sqrt(rgamma(n1, shape=shape, scale = scale))*sample(c(-1,1),n1,replace=T)
z2 = rnorm(n2)
zval = c(z1,z2)
pval = 2*pnorm(-abs(zval))


hpi0 = qvalue(pval)$pi0
lfdr.rst = mylocfdr(zval,pi0=hpi0)

path_lfdr = get_rej_path_lfdr(lfdr.rst$fdr)
path_qval = get_rej_path_pval(pval,hpi0)

plot(path_lfdr~path_qval)

pdf(file="lfdr_vs_qval_k_0_3.pdf", width=5, height=5.5,bg="white")
plot(path_lfdr ~ path_qval, pch=16,cex=0.5, xlab="rejection path of q-value", ylab = "rejection path of local fdr", main = "k = 0.3")
abline(0,1,lty=2, col="darkgreen")
dev.off()




# Hedenfalk data
library(qvalue)
data(hedenfalk)
pval = hedenfalk$p
length(pval)

zval = qnorm(pval/2)*sample(c(-1,1),length(pval),replace=T)
hpi0 = qvalue(pval)$pi0
lfdr.rst= mylocfdr(zval,hpi0)

path_lfdr = get_rej_path_lfdr(lfdr.rst$fdr)
path_qval = get_rej_path_pval(pval,hpi0)

plot(path_lfdr~path_qval)

pdf(file="lfdr_vs_qval_hedenfalk.pdf", width=5, height=5.5,bg="white")
plot(path_lfdr ~ path_qval, pch=16,cex=0.5, xlab="rejection path of q-value", ylab = "rejection path of local fdr")
abline(0,1,lty=2, col="darkgreen")
dev.off()







