library(qvalue)
library(ashr)
library(nleqslv)
library(locfdr)

z2_pdf<-function(y, sigma, pi){
  sum(sapply(seq(1:length(sigma)), function(x) pi[x]*dgamma(y,0.5,scale=2*(1+sigma[x]^2))))
}


z_pdf<-function(y, sigma, pi){
	  sum(sapply(seq(1:length(sigma)), function(x) pi[x]*dnorm(y,sd=sqrt(1+sigma[x]^2))))
}



z2_cdf<-function(y, sigma,pi){
	sum(sapply(seq(1:length(sigma)), function(x) pi[x]*pgamma(y,0.5,scale=2*(1+sigma[x]^2))))
}

z_cdf<-function(y,sigma,pi){
	sum(sapply(seq(1:length(sigma)), function(x) pi[x]*pnorm(y,sd=sqrt(1+sigma[x]^2))))
}


get_hat_pi0<-function(zval, lambda=0.5){
   thresh = qchisq(lambda,1)
   return(sum(zval^2<=thresh)/(0.5*length(zval)))

}


diag_quantile_z2_ash<-function(zval, ash.rst){
 
  library(nleqslv)

  N = length(zval)
  piv = ash.rst$fitted_g$pi
  sv = ash.rst$fitted_g$sd
  qv = seq(0.1,0.9,0.1)

  obs = quantile(zval^2, qv)

  quantile_eqn<-function(x){
    y = numeric(length(qv))
    for (i in 1:length(qv)){
       y[i] = z2_cdf(x[i],sv,piv) - qv[i]
    }
    y
  }


 sol =  nleqslv(obs, quantile_eqn)
 fitted = sol$x
 
 dv = sapply(fitted, function(x) z2_pdf(x, sv, piv))

 sdv = sqrt(qv*(1-qv)/(N*dv^2))

 tv = (obs - fitted)/sdv
 pvalue = 2*pnorm(-abs(tv))
 rst = rbind( obs, fitted, pvalue)
 if( length(which(pvalue<0.05 & fitted<obs)) >0){
	 cat(sprintf("\n *** Warning: potential violation of UA assumption *** \n\n"))
 }
 return(rst)
}



diag_quantile_z2_rbf<-function(zval, rbf.rst){
 
  library(nleqslv)

  N = length(zval)
  piv = rbf.rst$pi
  sv = rbf.rst$sd
  qv = seq(0.1,0.9,0.1)

  obs = quantile(zval^2, qv)

  quantile_eqn<-function(x){
    y = numeric(length(qv))
    for (i in 1:length(qv)){
       y[i] = z2_cdf(x[i],sv,piv) - qv[i]
    }
    y
  }


 sol =  nleqslv(obs, quantile_eqn)
 fitted = sol$x
 
 dv = sapply(fitted, function(x) z2_pdf(x, sv, piv))

 sdv = sqrt(qv*(1-qv)/(N*dv^2))

 tv = (obs - fitted)/sdv
 pvalue = 2*pnorm(-abs(tv))
 rbind(obs, fitted, pvalue)

}









diag_quantile_z_ash<-function(zval, ash.rst){
 
  library(nleqslv)
  N = length(zval)
  piv = ash.rst$fitted_g$pi
  sv = ash.rst$fitted_g$sd
  qv = seq(0.1,0.9,0.1)

  obs = quantile(zval, qv)

  quantile_eqn<-function(x){
    y = numeric(length(qv))
    for (i in 1:length(qv)){
       y[i] = z_cdf(x[i],sv,piv) - qv[i]
    }
    y
  }


 sol =  nleqslv(obs, quantile_eqn)
 fitted = sol$x
 
 dv = sapply(fitted, function(x) z_pdf(x, sv, piv))

 sdv = sqrt(qv*(1-qv)/(N*dv^2))

 tv = (obs - fitted)/sdv
 pvalue = 2*pnorm(-abs(tv))
 null = qnorm(qv,1)
 rbind(obs, fitted, pvalue)

}


diag_quantile_z_rbf<-function(zval, rbf.rst){
 
  library(nleqslv)
  N = length(zval)
  piv = rbf.rst$pi
  sv = rbf.rst$sd
  qv = seq(0.1,0.9,0.1)

  obs = quantile(zval, qv)

  quantile_eqn<-function(x){
    y = numeric(length(qv))
    for (i in 1:length(qv)){
       y[i] = z_cdf(x[i],sv,piv) - qv[i]
    }
    y
  }


 sol =  nleqslv(obs, quantile_eqn)
 fitted = sol$x
 
 dv = sapply(fitted, function(x) z_pdf(x, sv, piv))

 sdv = sqrt(qv*(1-qv)/(N*dv^2))

 tv = (obs - fitted)/sdv
 pvalue = 2*pnorm(-abs(tv))
 null = qnorm(qv,1)
 rbind(obs, fitted, pvalue)

}







### 06/23


diag_density_z2_call<-function(zval,sv, piv){

  zv2 = zval^2
  eds = density(zv2)
  max = max(eds$y)
 
  fds = sapply(zv2, function(x) z2_pdf(x,sv,piv))
  cutoff = min(zv2[which(fds<max)])

  x11()
  plot(fds[zv2>=cutoff]~zv2[zv2>=cutoff],col="brown",cex=0.25, ylab="density", xlab=expression(z^2))
  lines(eds,col="darkgreen")
  
}


## density plot
diag_density_z2_ash<-function(zval, ash.rst){

  piv = ash.rst$fitted_g$pi
  sv = ash.rst$fitted_g$sd
  
  diag_density_z2_call(zval, sv, piv)
  
}


diag_density_z2_rbf<-function(zval, rbf.rst){

  piv = rbf.rst$pi
  sv =  rbf.rst$sd
  
  diag_density_z2_call(zval, sv, piv)
  
}




diag_density_z_call<-function(zval, sv, piv){

  eds = density(zval)
  fds = sapply(zval, function(x) z_pdf(x,sv,piv))
  x11()
  plot(fds~zval,col="brown",cex=0.25, ylab="density", xlab=expression(z))
  lines(eds,col="darkgreen")

}


diag_density_z_ash<-function(zval, ash.rst){

  piv = ash.rst$fitted_g$pi
  sv = ash.rst$fitted_g$sd
  diag_density_z_call(zval, sv, piv)
}


diag_density_z_rbf<-function(zval, rbf.rst){

  piv = rbf.rst$pi
  sv =  rbf.rst$sd
  diag_density_z_call(zval, sv, piv)
}



 ### 06/25 

### check freq-Bayes FDR relationship



check_eqn_ash<-function(zval, ash.rst){
  N = length(zval)
  piv = ash.rst$fitted_g$pi
  sv = ash.rst$fitted_g$sd
  pip0 = get_lfdr(ash.rst)
  pval = sort(2*pnorm(-abs(zval)))
  fdp = cumsum(sort(pip0))/seq(1:N)
  nrej = max(which(fdp<=0.05))
  thresh = sort(zval^2,decreasing=T)[nrej]
  rst = c(nrej, pval[nrej]*get_pi0(ash.rst),sum(sort(pip0)[1:nrej])/N)
  names(rst) = c("rej", "freq", "bayes")
  return(rst)
}

diag_eqn_call <- function(zval, pip0, sv, piv){
  N = length(zval)	
  pval = sort(2*pnorm(-abs(zval)))
  yv = cumsum(sort(pip0))/N
  xv = piv[1]*pval
  x11()
  plot(yv ~xv, xlab = "frequentist FDR", ylab = "Bayesian FDR",  cex=0.2)
  abline(0,1, col="green", lty=2)
  dv = yv - xv
  x11()
  plot(dv ~ xv, xlab = "frequentist FDR", ylab = "BFDR - FDR", ylim = c(-max(abs(dv))-0.05, max(abs(dv))+0.05), cex= 0.2,col="brown")
  abline(h=0, col="green", lty=2)
}

diag_eqn <- function(zval, pip0, hat_pi0){
  N = length(zval)
  pval = sort(2*pnorm(-abs(zval)))
  yv = cumsum(sort(pip0))/N
  xv = hat_pi0*pval
  x11()
  plot(yv ~xv, xlab = "frequentist FDR", ylab = "Bayesian FDR",  cex=0.2)
  abline(0,1, col="green", lty=2)
  dv = yv - xv
  x11()
  plot(dv ~ xv, xlab = "frequentist FDR", ylab = "BFDR - FDR", ylim = c(-max(abs(dv))-0.05, max(abs(dv))+0.05), cex= 0.2,col="brown")
  abline(h=0, col="green", lty=2)
}


diag_eqn_truth<- function(zval, pip0, pi0){
  N = length(zval)
  pval = sort(2*pnorm(-abs(zval)))
  hpi0 = qvalue(pval)$pi0
  zv = hpi0*pval
  yv = cumsum(sort(pip0))/N
  xv = pi0*pval
  x11()
  plot(yv ~xv, xlab = "true FDR", ylab = "Bayesian FDR",  cex=0.2)
  points(zv ~xv, col="navyblue",  cex=0.2)
  abline(0,1, col="green", lty=2)
  dv = yv - xv
  x11()
  plot(dv ~ xv, xlab = "frequentist FDR", ylab = "BFDR - FDR", ylim = c(-max(abs(dv))-0.05, max(abs(dv))+0.05), cex= 0.2,col="brown")
  abline(h=0, col="green", lty=2)
}



diag_eqn_ash<-function(zval, ash.rst){
  piv = ash.rst$fitted_g$pi
  sv = ash.rst$fitted_g$sd
  pip0 = get_lfdr(ash.rst)
  diag_eqn_call(zval, pip0,  sv, piv)
}

diag_eqn_rbf<-function(zval, rbf.rst){
  piv = rbf.rst$pi
  sv =  rbf.rst$sd
  pip0 = rbf.rst$pip0
  diag_eqn_call(zval, pip0,  sv, piv)
}

diag_eqn_locfdr<-function(zval, locfdr.rst){
   pi0_est = locfdr.rst$fp0[1,3]
   piv = c(pi0_est)
   diag_eqn_call(zval, locfdr.rst$fdr, c(0), piv)
}


fdr_control<-function(pip0, alpha){

 rst = list( critical_val = thresh, rej_quantile = round(100*(1-rej_bayes/length(pip0)),1), rej_set = index)
  return(rst)


}


check_truth_pip0<-function(pip0, truth, alpha = 0.05){
  
  pip0_sort = sort(pip0)
  bfdp = cumsum(pip0_sort)/seq(1:length(pip0))
  rej_bayes = max(which(bfdp<=alpha))
  thresh = pip0_sort[rej_bayes]
  index = which(pip0<= thresh)
  
  cat(sprintf("\nrejection quantile: >= %s%s\n",round(100*(1-rej_bayes/length(pip0)),1),"%" ))
  cat(sprintf("critical value (lfdr): %s\n\n", round(thresh,3)))

  if(length(index) == 0){
	fdr = 0
	power = 0
  }else{
       	false_rej = length(which(truth[index]==0)) 
        fdr = false_rej/length(index)
	power = length( which(truth[index] == 1))/sum(truth)
  }
  rst = c(fdr, power)
  names(rst) = c("fdr", "power")
  return(rst)
}

check_truth_qval<-function(zval, truth, alpha=0.05){
  pval = 2*pnorm(-abs(zval))
  qrst = qvalue(pval,fdr.level=alpha)

  index = which(qrst$sig)

  cat(sprintf("\nrejection quantile: >= %s%s\n",round(100*(1-length(index)/length(pval)),1),"%" ))
  cat(sprintf("critical value (pvalue): %e\n\n",sort(pval)[length(index)] ))
  



  if(length(index) == 0){
     fdr = 0
     power = 0
   }else{
     false_rej = length(which(truth[index]==0))
     fdr = false_rej/length(index)
     power = length( which(truth[index] == 1))/sum(truth)
   }
   rst = c(fdr, power)
   names(rst) = c("fdr", "power")
   return(rst)
}  



