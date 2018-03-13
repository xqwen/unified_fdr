compute_log10_BF<-function(zval, k){
      return(-0.5*log10(1+k)+0.5*(k/(1+k))*zval^2/log(10))
}

compute_log10_lik<-function(pi){
	return(sum(log10((pi+(1-pi)*BFv))))
}

compute_pip0<-function(pi){
	return( pi/(pi+(1-pi)*BFv))
}

pi0 = 0.5 
K = 10


set.seed(123)

N = 200
n1 = round((1-pi0)*N)
n2 = N-n1
z1 = rnorm(n1,sd=sqrt(1+K))
z2 = rnorm(n2)
zval = c(z1,z2)
BFv = 10^compute_log10_BF(zval,K)
pip0 = compute_pip0(pi0)

bfdr = cumsum(sort(pip0))/seq(1:N)

pval = 2*pnorm(-abs(zval))
ffdr = pi0*sort(pval)*N/seq(1:N)
pdf(file="bfdr_vs_fdr_m_200.pdf", pointsize = 16, width=7, height=7.5, bg="white")
plot(bfdr~ffdr,pch=16, cex=0.7, xlim=c(0,0.53), ylim = c(0,0.53), xlab="FDR", ylab = "BFDR", main = "m=200")
abline(0,1, col="darkgreen")
dev.off()
sqrt(mean((ffdr-bfdr)^2))




N = 2000
n1 = round((1-pi0)*N)
n2 = N-n1
z1 = rnorm(n1,sd=sqrt(1+K))
z2 = rnorm(n2)
zval = c(z1,z2)
BFv = 10^compute_log10_BF(zval,K)
pip0 = compute_pip0(pi0)

bfdr = cumsum(sort(pip0))/seq(1:N)

pval = 2*pnorm(-abs(zval))
ffdr = pi0*sort(pval)*N/seq(1:N)

pdf(file="bfdr_vs_fdr_m_2k.pdf", pointsize = 16, width=7, height=7.5, bg="white")
plot(bfdr~ffdr,pch=16, cex=0.7, xlab="FDR", ylab = "BFDR", main = "m=2000")
abline(0,1, col="darkgreen")
dev.off()
sqrt(mean((ffdr-bfdr)^2))


N = 20000
n1 = round((1-pi0)*N)
n2 = N-n1
z1 = rnorm(n1,sd=sqrt(1+K))
z2 = rnorm(n2)
zval = c(z1,z2)
BFv = 10^compute_log10_BF(zval,K)
pip0 = compute_pip0(pi0)

bfdr = cumsum(sort(pip0))/seq(1:N)

pval = 2*pnorm(-abs(zval))
ffdr = pi0*sort(pval)*N/seq(1:N)
pdf(file="bfdr_vs_fdr_m_20k.pdf", pointsize = 16, width=7, height=7.5, bg="white")
plot(bfdr~ffdr,pch=16, cex=0.7, xlab="FDR", ylab = "BFDR", main = "m=20000")
abline(0,1, col="darkgreen")
dev.off()
sqrt(mean((ffdr-bfdr)^2))




