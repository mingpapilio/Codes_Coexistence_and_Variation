library(gtools)

## The thermal perfornamce curve #####
  x_opt1<- 30
  x_opt2<- 17
  x_max1<- 35
  x_max2<- 23
  sigma1<- 5
  sigma2<- 2
  pcurv_1<- function(x){
    if(x<x_opt1) return(exp(-((x-x_opt1)/2/sigma1)^2))
    if(x_max1>x && x>=x_opt1) return(1-((x-x_opt1)/(x_opt1-x_max1))^2)
    if(x>=x_max1) return(0)
  }
  pcurv_2<- function(x){
    if(x<x_opt2) return(exp(-((x-x_opt2)/2/sigma2)^2))
    if(x_max2>x && x>=x_opt2) return(1-((x-x_opt2)/(x_opt2-x_max2))^2)
    if(x>=x_max2) return(0)
  }
  

# Read in the files #####

raw<-read.csv("out.txt",sep="\t")

dev.new()
par(mar=c(5,5,2,2))
plot(c(0,0), xlim=c(1, 1000), ylim=c(-10, 50), col="white", xlab="window", ylab="averages")
ACF<- acf(raw$env, lag.max= 1000,plot=F)

for(i in 1: 100){
  span= 10*i
  tmp_y<- sapply(split(raw$env, (seq_along(raw$env) - 1) %/% span),mean)
  tmp_x<- rep(span, length(tmp_y))
  tmp<- data.frame()
  tmp<- cbind(tmp_x, tmp_y)
  points(tmp, pch=20)
  avg<- mean(tmp_y)
  sd<- sd(tmp_y)
  points(x=span,y=avg+2*sd, pch="-", col="red", cex=2)
  points(x=span,y=avg-2*sd, pch="-", col="red", cex=2)
  tmp_1= tmp_2= 0
  for(j in 1:length(tmp_y)){
    tmp_1= tmp_1+ pcurv_1(tmp_y[j])/length(tmp_y)
    tmp_2= tmp_2+ pcurv_2(tmp_y[j])/length(tmp_y)
  }
  points(x=span,y=tmp_1*50*ACF$acf[span+1], pch="-", col="orange", cex=2)
  points(x=span,y=tmp_2*50*ACF$acf[span+1], pch="-", col="skyblue", cex=2)
}

aa<- fft(raw$env)
aa.spec<- spectrum(aa,log="no",plot=F)
dev.new()
par(mar=c(5,5,2,2))
plot(aa.spec$spec~aa.spec$freq,xlab="frequency",ylab="spectral density",log="xy",type="l")

dev.new()
par(mar=c(5,5,2,2))
plot(spectrum(raw$env[1:1000],log="no",plot=F), ylim=c(5E-5, 5E4),
     ylab="amplitude", main="")
# Time series #####
raw<-read.csv("out.txt",sep="\t")
dev.new()
par(mar=c(5,5,2,2))
plot(raw$env[1:1000], ylab="temperature",xlab="time unit", ylim=c(-10,50),type="l")

#####
## Start read in
raw<-read.csv("out.txt",sep="\t")
# raw$env[1]= 0
# Original time series #####
ACF<- acf(raw$env, lag.max= 70,plot=F)
ACF<- ACF[1:70]
dev.new()
par(mar=c(5,5,2,2))
plot(ACF, main="", ylim=c(0,1))
PACF<- pacf(raw$env, lag.max= 71,plot=F)
dev.new()
par(mar=c(5,5,2,2))
plot(PACF, main="", ylim=c(0,1))

#####
raw<-read.csv("out.txt",sep="\t")
ACF<- acf(raw$env, lag.max= 5,plot=F)
span<- 1
ACF<- ACF$acf[span+1]
tmp_y<- sapply(split(raw$env, (seq_along(raw$env) - 1) %/% span),mean)
for(j in 1:length(tmp_y)){
  tmp_1= tmp_1+ pcurv_1(tmp_y[j])/length(tmp_y)
  tmp_2= tmp_2+ pcurv_2(tmp_y[j])/length(tmp_y)
}
c(tmp_1*ACF, tmp_2*ACF)

