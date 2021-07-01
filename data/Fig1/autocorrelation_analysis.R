
## Reference time series: white noise
wn<- rnorm(30001, 0, 15)
acf(wn, lag.max=30, plot=T)
sine<- sin(seq(0,10*pi,length.out = 3000))
  dev.new()
  par(mfrow=c(2,1))
  ACF<- acf(sine, lag.max= 25,plot=F)
  PACF<- pacf(sine, lag.max= 25,plot=F)
  plot(ACF, main="Time series of environmantal conditions")
  plot(PACF, main="Time series of environmantal conditions")

## Start read in
raw<-read.csv("out.txt",sep="\t")
# raw$env[1]= 0
# Original time series #####
  dev.new()
  par(mfrow=c(2,1))
  ACF<- acf(raw$env, lag.max= 10,plot=F)
  ACF<- ACF[1:10]
  PACF<- pacf(raw$env, lag.max= 10,plot=F)
  plot(ACF, main="Time series of environmantal conditions")
  plot(PACF, main="Time series of environmantal conditions")
# First-order derivative #####
env.diff<- diff(raw$env)
  dev.new()
  par(mfrow=c(2,1))
  ACF<- acf(env.diff, lag.max= 100)
  PACF<- pacf(env.diff, lag.max= 100)
  plot(ACF, main="Time series of first order derivative")
  plot(PACF, main="Time series of first order derivative")

# Second-order derivative #####
env.diff.2<- diff(diff(raw$env))
  dev.new()
  par(mfrow=c(2,1))
  ACF<- acf(env.diff.2, lag.max= 100)
  PACF<- pacf(env.diff.2, lag.max= 100)
  plot(ACF, main="Time series of second order derivative")
  plot(PACF, main="Time series of second order derivative")
  
### Estimating the AR model by Yule-Walker method #####
env.yw<- ar.yw(raw$env, partialacf= T
               , order.max= 72
              )
  ### mean estimate
  env.yw$x.mean
  ### phi estimates
  env.yw$ar
  ### their SD
  sqrt(diag(env.yw$asy.var.coef))
  ### error variance estimate
  env.yw$var.pred
###
dev.new()
hist(env.yw$resid)
qqnorm(env.yw$resid)
shapiro.test(env.yw$resid)

### 
env.mle<- ar.mle(raw$env, order.max= 72, partialacf= T)

###
env.fit1<- arima(raw$env, order= c(2,0,0), seasonal= list(order=c(2,0,0), period= 71), method= "CSS")
tsdiag(env.fit1, gof.lag= 280)

env.fit2<- arima(raw$env, order= c(2,1,1), seasonal= list(order=c(2,1,1), period= 70), method= "CSS")
tsdiag(env.fit2, gof.lag= 280)

env.fit3<- arima(raw$env, order= c(71,1,1), method= "CSS")
tsdiag(env.fit3, gof.lag= 280)
  
### Prediction
env.pr<- predict(env.yw, n.ahead= 150)
pr.upr= env.pr$pred+ env.pr$se
pr.lwr= env.pr$pred- env.pr$se
time= 15001: 20151
dev.new()
plot(time, raw$env[time], type="l", xlim= c(15000, 20200))
  lines(env.pr$pred, col="red", type="l")
  lines(pr.upr, col="grey", lty="dashed")
  lines(pr.lwr, col="grey", lty="dashed")

#####
  x<- ts(cbind(raw$N1, raw$N2, raw$env))
  s<- spec.pgram(x)
  
  x<- ts(cbind(raw$N1, raw$env))
  s<- spec.pgram(x, kernel("modified.daniell",c(5,5)), log="no")
  f<- qf(.999, 2, s$df-2)
  c<- f/(21+f)
  dev.new()
  plot(s, plot.type="coh", ci.lty=2)
  abline(h=c)
  abline(v=1/70)

k<- kernel("modified.daniell", 6)
envf<- kernapply(raw$env, k)
spectrum(envf, spans=9, log="no") # Has to include "no" log
