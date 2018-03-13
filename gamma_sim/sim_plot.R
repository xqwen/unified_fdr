d = read.table("sim.summary")
attach(d)

pdf(file="est_pi0_za_vs_ua.pdf", height=6,width=8,bg="white")
V4[V4>1]=1
V5[V5>1]=1
boxplot(V4~V2, ylim=c(0,1), col="navyblue",outline=FALSE,ylab=expression(paste(pi[0]," estimate")),xlab = "shape parameter")
boxplot(V5~V2, ylim=c(0,1), col="darkgreen",add=T,outline=FALSE)
abline(h=0.6,lty=2)
legend("bottomleft", c("UA estimates", "ZA estimates"), fill =c("darkgreen","navyblue"),cex=1)
dev.off()


