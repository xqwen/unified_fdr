library(qvalue)
library(ashr)
library(nleqslv)



############# functions for posterior predictive checking #########################


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



model_check_ash<-function(zval, ash.rst){

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
    status = 0
    if( length(which(pvalue<0.05 & fitted<obs)) >0){
        cat(sprintf("\n *** Warning: potential violation of UA assumption *** \n\n"))
        status=1
    }
    return(list(rst,status))
}



model_check_rbf<-function(zval, rbf.rst){

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




####################### utility functions ########################

get_hat_pi0<-function(zval, lambda=0.5){
    thresh = qchisq(lambda,1)
    return(sum(zval^2<=thresh)/(0.5*length(zval)))

}



##################### rejection path checking #####################


get_rej_path_lfdr<-function(pip0_vec){
    return(cumsum(sort(pip0_vec))/(1:length(pip0_vec)))
}

get_rej_path_pval<-function(pval_vec, hat_pi0){
    return(sort(pval_vec)*hat_pi0*length(pval_vec)/(1:length(pval_vec)))
}


plot_rejection_path <- function(path1, path2, name = c("Rejection path 1", "Rejection path 2")){
    rp.plot = plot(path1 ~ path2, xlab = name[2], ylab = name[1],  cex=0.2)
    abline(0,1, col="green", lty=2)
}






