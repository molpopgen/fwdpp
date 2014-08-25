library(ggplot2)
library(dplyr)

thetas=array()
times=array()
mems=array()
version=array()
I=1
NREPS=128
for( theta in c(10,50,100,250,500))
    {
        for(i in 1:NREPS)
            {
                fn=paste("time.",i,".txt",sep="")
                fn.cpp=paste("cpp11.",theta,"/",fn,sep="")
                if(file.exists(fn.cpp))
                    {
                        x=scan(fn.cpp,quiet=T)
                        if(length(x)>0)
                            {
                                thetas[I]=theta
                                version[I]="0.2.5"
                                times[I]=as.numeric(x[1])
                                mems[I]=as.numeric(x[2])
                            }
                    }
                I=I+1
            }
        for(i in 1:NREPS)
            {
                fn=paste("time.",i,".txt",sep="")
                pub.cpp=paste("pub.",theta,"/",fn,sep="")
                if(file.exists(pub.cpp))
                    {
                        x=scan(pub.cpp,quiet=T)
                        if(length(x)>0)
                            {
                                thetas[I]=theta
                                version[I]="0.2.4"
                                times[I]=as.numeric(x[1])
                                mems[I]=as.numeric(x[2])
                            }
                    }
                I=I+1
            }
    }

d=data.frame(cbind(thetas,version,times,mems))

write.table(d,file="time_mem.txt",quote=F,row.names=F,col.names=T)
d=read.table("time_mem.txt",header=T)

d.groups = group_by(d,thetas,version)
d.summ = summarise(d.groups,mm=mean(mems),mt=mean(times),mdt=median(times),mint=min(times),maxt=max(times))
#d.summ = ddply(d,c("thetas","version"),summarise,mt=mean(as.numeric(times)),mm=mean(as.numeric(mems)))
print(d.summ)

d.summ.plot=ggplot(d.summ,aes(x=thetas,y=mt/60)) + geom_point(aes(color=version,shape=version)) + geom_line(aes(color=version)) + xlab(expression(paste(theta," = ",rho))) + ylab("Mean run time (minutes)")

#d.summ.plot + scale_y_continuous(trans="log")

d.summ.mem.plot=ggplot(d.summ,aes(x=thetas,y=mm/(4*1024))) + geom_point(aes(color=version,shape=version)) + geom_line(aes(color=version)) + xlab(expression(paste(theta," = ",rho))) + ylab("Mean peak RAM use (megabytes)")
#d.summ.mem.plot + scale_y_continuous(trans="log")
