## Time series analysis ##

## Read in the files #####
raw<-read.csv("out.txt",sep="\t")

## Plot temperature time series #####
  dev.new()
  par(mar=c(5,5,2,2))
  plot(raw$env[1:1000], ylab="temperature",xlab="time unit", ylim=c(-10,50),type="l")
  
## Plot the population dynamics #####
  dev.new()
  par(mar=c(5,5,2,2))
  plot(raw$N1[1:1000], ylab="population size",xlab="time unit", ylim=c(0,10000),type="l",col="orange")
  points(raw$N2[1:1000], type="l", col="skyblue")
  
## Plot the windowed averages #####
  dev.new()
  par(mar=c(5,5,2,2))
  plot(c(0,0), xlim=c(1, 200), ylim=c(-10, 50), col="white", xlab="window", ylab="averages")
  
  for(i in 1: 50){
    span= 4*i
    tmp_y<- sapply(split(raw$env, (seq_along(raw$env) - 1) %/% span),mean)
    tmp_x<- rep(span, length(tmp_y))
    tmp<- data.frame()
    tmp<- cbind(tmp_x, tmp_y)
    points(tmp, pch=20)
    avg<- mean(tmp_y)
    sd<- sd(tmp_y)
    points(x=span,y=avg+2*sd, pch="-", col="red", cex=2)
    points(x=span,y=avg-2*sd, pch="-", col="red", cex=2)
  }
## Autocorrelation analysis #####
  ACF<- acf(raw$env, lag.max= 70,plot=F)
  ACF<- ACF[1:70]
  dev.new()
  par(mar=c(5,5,2,2))
  plot(ACF, main="", ylim=c(0,1))
  # PACF
  PACF<- pacf(raw$env, lag.max= 71,plot=F)
  dev.new()
  par(mar=c(5,5,2,2))
  plot(PACF, main="", ylim=c(0,1))
  
## Fast Fourier transformation #####
  dev.new()
  par(mar=c(5,5,2,2))
  plot(spectrum(raw$env[1:1000],log="no",plot=F), ylim=c(5E-5, 5E4),
       ylab="amplitude", main="")

## The thermal perfornamce curve #####
  x_opt1<- 30
  x_opt2<- 17
  x_max1<- 35
  x_max2<- 23
  sigma1<- 5
  sigma2<- 2
  dev.new()
  par(mar=c(5,5,2,2))
  curve((x<x_opt1)*(exp(-((x-x_opt1)/2/sigma1)^2))
        +(x_max1>x & x>=x_opt1)*(1-((x-x_opt1)/(x_opt1-x_max1))^2)
        +(x>=x_max1)*(0), lwd=2, n=200,
        xlim=c(-10,50),ylim=c(0,1),col="orange",
        ylab="Thermal performance function")
  curve((x<x_opt2)*(exp(-((x-x_opt2)/2/sigma2)^2))
        +(x_max2>x & x>=x_opt2)*(1-((x-x_opt2)/(x_opt2-x_max2))^2)
        +(x>=x_max2)*(0)
        ,col="skyblue", lwd=2, n=200, add=T)
  
##