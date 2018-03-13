library(qvalue)


mylocfdr<-function(zval,pi0 = -1, truncation = TRUE){
    if(pi0 == -1){	
       pi0 = qvalue(2*pnorm(-abs(zval)))$pi0
    }
     df = approxfun(density(zval))
     pip0 = pi0*dnorm(zval)/df(zval)
     if(truncation == TRUE){ 
	     pip0[pip0>1] = 1
     }
     return(list(pi0 = pi0, fdr=pip0))
}   


#convert_lfdr<-function(pval, locfdr.rst, hpi0 = -1){
#   pip0 = locfdr.rst$fdr
#   if(hpi0 == -1){
#        hpi0 =  pi0_qval = qvalue(pval)$pi0
#   }
#   pi0_est = locfdr.rst$fp0[1,3]
#   return(locfdr.rst$fdr*hpi0/pi0_est)
#}

