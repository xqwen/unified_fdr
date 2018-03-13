source("mylocfdr.R")

sim_data<-function(shape){
    N = 20000
    pi0 = 0.55
    n1 = round((1-pi0)*N)
    n2 = N-n1

    scale = 10
    #set.seed(123)

    z1 = sqrt(rgamma(n1, shape=shape, scale = scale))*sample(c(-1,1),n1,replace=T)
    z2 = rnorm(n2)
    zval = c(z1,z2)
    pval = 2*pnorm(-abs(zval))


    lfdr.rst = mylocfdr(zval)

    cor(lfdr.rst$fdr, pval, method="spearman")
}



set.seed(123)
sv = seq(0.05, 1.00,0.01)
simv = sapply((1:50), function(x) sapply(sv,sim_data))
mv = apply(simv,1,mean)
plot(mv ~ sv, pch=16, xlab="shape parameter (k)", ylab = "rank correlation")


