data<-read.csv("out.txt",sep="\t",skip=2)
cor(data[,3],data[,4])

data<-read.csv("out.txt",sep="\t",skip=2)
aa<- fft(data[,4])
aa.spec<- spectrum(aa,log="no",plot=F)
dev.new()
plot(aa.spec$spec~aa.spec$freq,xlab="frequency",ylab="spectral density",log="xy",type="l")


data<-read.csv("MpalaStation_tmp.txt",sep="\t",skip=2)
aa<- fft(data[,1])
# dev.new()
# spectrum(aa,log="yes")
aa.spec<- spectrum(aa,log="no",plot=F)
dev.new()
plot(aa.spec$spec~aa.spec$freq,xlab="frequency",ylab="spectral density",log="xy",type="l")


data<-read.csv("Toolik_tmp.txt",sep="\t",skip=2)
xx<- c(data[1:40,1])
aa<- fft(xx)
aa<- fft(data[,1])
# dev.new()
# spectrum(aa,log="yes")
aa.spec<- spectrum(aa,log="no",plot=F)
dev.new()
plot(aa.spec$spec~aa.spec$freq,xlab="frequency",ylab="spectral density",log="xy",type="l")


t <- seq(0,1024,by=0.1)
x <- cos(2*pi*t) + 0.75*sin(2*pi*4*t) + 2*sin(2*pi*6*t)
spectrum(x,log="no")#,span=5,plot=FALSE)

