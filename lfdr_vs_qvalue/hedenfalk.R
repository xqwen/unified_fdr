source("mylocfdr.R")

convert_lfdr<-function(lfdr_vec,pi0 = -1){
    if(pi0<=0 || pi0 >=1){
       pi0 = mean(lfdr_vec)
    }
    minp = min(lfdr_vec[lfdr_vec>0])
    lfdr_vec[lfdr_vec == 0] = 0.1*minp

    lfdr_vec_sort = sort(lfdr_vec)
    bfdr = cumsum(lfdr_vec_sort)/(1:length(lfdr_vec_sort))
    pval = (1:length(lfdr_vec))*bfdr/(length(lfdr_vec)*pi0)
    pv = rep(0,length(lfdr_vec))
    pv[order(lfdr_vec)] = pval
    return(2)
}




data(hedenfalk)
pval = hedenfalk$p
length(pval)

zval = qnorm(pval/2)*sample(c(-1,1),length(pval),replace=T)
lfdr.rst= mylocfdr(zval)

pval2 = convert_lfdr(lfdr.rst$fdr)
