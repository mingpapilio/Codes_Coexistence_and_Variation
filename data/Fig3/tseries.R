## Plot the time series
raw<-read.csv("out.txt",sep="\t")
dev.new()
par(mar=c(5,4,2,2))
plot(raw[,3],raw[,1],type="l",lwd=2,col="orange",xlab="Time",ylab="Population size",
     ylim=c(0,10000))
points(raw[,3],raw[,2],type="l",lwd=2,col="skyblue")