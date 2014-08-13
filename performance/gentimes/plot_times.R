x1=scan("times")
x2=scan("times2")
x3=scan("times3")

png("gentimes.png")
plot(ksmooth(1:length(x3),x3,bandwidth=100),type="l",xlab="Generation",ylab="Run time per gen (seconds, smoothed in 100 generation intervals)")
lines(ksmooth(1:length(x2),x2,bandwidth=100),type="l",col="blue")
lines(ksmooth(1:length(x1),x1,bandwidth=100),type="l",col="red")
dev.off()
