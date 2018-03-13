library(qvalue)
data(hedenfalk)
pval = hedenfalk$p
m = length(pval)
m
qval_rst = qvalue(pval,fdr.level=0.05)
pi0 = qval_rst$pi0
path_qval = pi0*sort(pval)*m/seq(1:m)
max(which(path_qval<0.05))
path_BH = sort(pval)*m/seq(1:m)
pdf(file="rej_path.pdf", width=5, height=5.5,bg="white")
plot(path_BH~path_qval, xlim=c(0,1), ylim=c(0,1), pch=16,cex=0.5, xlab="rejection path of q-value", ylab = "rejection path of B-H")
abline(0,1, lty=2, col="green")
dev.off()
