# usage: Rscript sim.R shape index



source("ufdr.R")

args = commandArgs(trailingOnly = TRUE)

pi0 = 0.60
shape = as.numeric(args[1])
index = as.numeric(args[2])

set.seed(index)

N = 10000
n1 = round((1-pi0)*N)
n2 = N-n1

scale  = runif(1,min=5, max=15)
#z1 = rnorm(n1, sd=sqrt(1+K))
z1 = sqrt(rgamma(n1, shape=shape, scale = scale))
z2 = rnorm(n2)

zval = c(z1,z2)
ash.rst = ash(zval, rep(1,length(zval)),mixcomp="normal")
ash_pi0 = get_pi0(ash.rst)
hat_pi0 = get_hat_pi0(zval)
ash_diag = model_check_ash(zval, ash.rst)
status = ash_diag[[2]]
outd = c(round(pi0,2), round(shape,3), round(scale,3), round(hat_pi0,3), round(ash_pi0,3), round(status,0))
write(file=paste("sim_out/sim.", index,".rst",sep=""), outd, ncol=7)




