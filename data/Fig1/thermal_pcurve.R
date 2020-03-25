## Plot the thermal perfornamce curve
  x_opt1<- 30
  x_opt2<- 17
  x_max1<- 35
  x_max2<- 23
  sigma1<- 5
  sigma2<- 2
  dev.new()
  curve((x<x_opt1)*(exp(-((x-x_opt1)/2/sigma1)^2))
        +(x_max1>x & x>=x_opt1)*(1-((x-x_opt1)/(x_opt1-x_max1))^2)
        +(x>=x_max1)*(0)
        , xlim=c(0,40),ylim=c(0,1),col="orange"
        , ylab="Thermal performance function")
  curve((x<x_opt2)*(exp(-((x-x_opt2)/2/sigma2)^2))
        +(x_max2>x & x>=x_opt2)*(1-((x-x_opt2)/(x_opt2-x_max2))^2)
        +(x>=x_max2)*(0)
        ,col="skyblue",add=T)
