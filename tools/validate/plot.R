n=commandArgs(trailing=TRUE)
f=read.table(n[1],header=T)
m=read.table(n[2],header=T)

pdf(n[3],pointsize=20,height=10,width=10)
par(mfrow=c(3,2))

plot(ecdf(f$S),main="S",do.points=F)
plot(ecdf(m$S),add=TRUE,col="red",do.points=F)

plot(ecdf(f$pi),main="pi",do.points=F)
plot(ecdf(m$pi),add=TRUE,col="red",do.points=F)


plot(ecdf(f$next.),main="derived singletons",do.points=F)
plot(ecdf(m$next.),add=TRUE,col="red",do.points=F)

plot(ecdf(f$hdiv),main="hdiv",do.points=F)
plot(ecdf(m$hdiv),add=TRUE,col="red",do.points=F)

plot(ecdf(f$nhaps),main="# Haplotypes",do.points=F)
plot(ecdf(m$nhaps),add=TRUE,col="red",do.points=F)

plot(ecdf(f$rm),main="minrec (H&K)",do.points=F)
plot(ecdf(m$rm),add=TRUE,col="red",do.points=F)
dev.off()

##p-values
s.pv = ks.test(f$S,m$S)$p.value
pi.pv = ks.test(f$pi,m$pi)$p.value
next.pv = ks.test(f$next.,m$next.)$p.value
hdiv.pv = ks.test(f$hdiv,m$hdiv)$p.value
nhaps.pv = ks.test(f$nhaps,m$nhaps)$p.value
rm.pv = ks.test(f$rm,m$rm)$p.value

print(paste(s.pv,pi.pv,next.pv,hdiv.pv,nhaps.pv,rm.pv,sep=" "))
