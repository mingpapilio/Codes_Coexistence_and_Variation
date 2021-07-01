library(plotly)
library(viridis)
## List of variables =====
  list_mean<- seq(2, 38, by=3)
  list_long<- seq(0.5, 20.5, by=5)
  list_short<- seq(20.5, 0.5, by= -5)
  table_list_mean<- c(rep(2, 5), rep(5,5), rep(8,5), rep(11,5), rep(14,5), rep(17,5), rep(20,5), rep(23,5), rep(26,5), rep(29,5), rep(32,5), rep(35,5), rep(38,5))


## Plotting =====
rr<- rep(NA, 11)
gg= bb= rr
data<-read.csv("summary.txt",sep="\t")
rep<- 100

rr[1:6]<- round(seq(86, 202, length.out= 6))
gg[1:6]<- round(seq(53, 122, length.out= 6))
bb[1:6]<- round(seq(46, 44, length.out= 6))

rr[6:11]<- round(seq(202, 250, length.out= 6))
gg[6:11]<- round(seq(122, 214, length.out= 6))
bb[6:11]<- round(seq(44, 137, length.out= 6))

trace1<- list(
  z= matrix(data$coexist,nrow=5,ncol=13),
  type= "contour"
  ,colorscale= cbind(seq(0,1,by=1/10), viridis(11))
  ## viridis(), inferno() are recommended
)
p <- plot_ly(
  contours = list(
    start= 0.1*rep,
    end= 0.9*rep,
    size= 0.1*rep
  ))
p<- add_trace(p, z=trace1$z,x=list_mean,y=list_long, 
              colorscale=trace1$colorscale, type=trace1$type,
              colorbar=list(tickfont=list(size=18,color='black')))
p<- layout(p, xaxis=list(tickfont=list(size=24)),yaxis=list(tickfont=list(size=24)))
p
