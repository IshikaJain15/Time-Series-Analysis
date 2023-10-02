#ISHIKA JAIN
#r0915387
#Advanced Time Series Analysis Homework


#install packages
install.packages("rgl")
install.packages("MASS")
install.packages("CADF")
install.packages("fGarch")
install.packages("forecast")
install.packages("fpp")

#
rm(list=ls())
#load data
mydata <- read.table(file=file.choose(),header=TRUE)
attach(mydata)
exts <- ts(EXINUS, frequency = 12, start=c(2000,1))
dexts <- diff(exts)
logexts<-ts(log(EXINUS), frequency = 12)
dlogexts<-diff(logexts)
#libraries

library(CADFtest)
library(fGarch)
library(forecast)
library(vars)
library(fpp)
library(rgl)

#plot data

plot.ts(exts)
plot.ts(dexts)
plot.ts(logexts)
plot.ts(dlogexts)
par(mfrow=c(2,2))
#correlograms
acf(exts)
acf(dexts) #MA(1)
acf(logexts)
acf(dlogexts)#ma1
pacf(exts)
pacf(dexts)
pacf(logexts)
pacf(dlogexts)
plot(dexts)
monthplot(dexts)
monthplot(dlogexts)
#Tests for stationarity and white noise
CADFtest(exts, type= "drift", criterion= "BIC", max.lag.y=max.lag)  #do not reject, not stationary
Box.test(exts, lag = max.lag, type = "Ljung-Box") #reject, not white noise
max.lag<-round(sqrt(length(exts)))
CADFtest(dexts, type= "drift", criterion= "BIC", max.lag.y=max.lag) #reject CADFtest, stationary
Box.test(dexts, lag = max.lag, type = "Ljung-Box") #reject, not white noise

CADFtest(dlogexts, type= "drift", criterion= "BIC", max.lag.y=max.lag) #reject CADFtest, stationary
Box.test(dlogexts, lag = max.lag, type = "Ljung-Box") #reject, not white noise
par(mfrow=c(2,2))

#Intrayear pattern

one.year <- rep(1:12, 22) 
oneyear <- append(one.year, c(1, 2, 3, 4, 5, 6, 7, 8, 9))
D <- 12; sp <- 2*pi*oneyear/D
intrayear.model <- lm(exts ~ sin(sp)+cos(sp)+sin(2*sp)+cos(2*sp))
summary(intrayear.model)


#Estimate MA(1) ARIMA(0,1,1)
par(mfrow=c(1,2))
arima1 <- arima(exts, order = c(0,1,1))
arima1 #significant coefficient
acf(arima1$res, main="Residuals- ARIMA(0,1,1)")
acf(arima1$res^2, main="Squared Residuals ARIMA(0,1,1)") #heteroskedasticity might be present
Box.test(arima1$res, lag = max.lag, type = "Ljung-Box") #do not reject, white noise
fcstarima1 <- predict(arima1, n.ahead = 12)
expected <- fcstarima1$pred
lower <- fcstarima1$pred-qnorm(0.975)*fcstarima1$se
upper <- fcstarima1$pred+qnorm(0.975)*fcstarima1$se
plot.ts(exts, xlim=c(2021,2024), ylim=c(60,90), main="ARIMA(0,1,1)")
lines(expected,col="red")
lines(lower,col="blue")
lines(upper,col="blue")



#Estimate ARIMA(1,0,1)

arima2 <- arima(exts, order = c(1,0,1))
arima2 #both terms significant
plot.ts(arima2$res)
acf(arima2$res, main="ARIMA(1,0,1)")
Box.test(arima2$res, lag = max.lag, type = "Ljung-Box") #do not reject, white noise
fcstarima2 <- predict(arima2, n.ahead = 12)
expected2 <- fcstarima2$pred
lower2 <- fcstarima2$pred-qnorm(0.975)*fcstarima2$se
upper2 <- fcstarima2$pred+qnorm(0.975)*fcstarima2$se
plot.ts(exts,xlim=c(2021, 2024),ylim=c(60,90), main="ARIMA(1,0,1)")
lines(expected2,col="red")
lines(lower2,col="blue")
lines(upper2,col="blue")

#Estimate ARIMA(1,1,0)

arima3 <- arima(exts, order = c(1,1,0))
arima3 # AR term significant
plot.ts(arima3$res)
acf(arima3$res, main="AR(1) in differences- ARIMA(1,1,0)")
Box.test(arima3$res, lag = max.lag, type = "Ljung-Box") #do not reject, white noise
fcstarima3 <- predict(arima3, n.ahead = 12)
expected3 <- fcstarima3$pred
lower3 <- fcstarima3$pred-qnorm(0.975)*fcstarima3$se
upper3 <- fcstarima3$pred+qnorm(0.975)*fcstarima3$se
plot.ts(exts,xlim=c(2021, 2024),ylim=c(60,90), main="ARIMA(1,1,0)")
lines(expected3,col="red")
lines(lower3,col="blue")
lines(upper3,col="blue")



#GARCH effects

mymodel<-garchFit(~arma(0,1) + garch(1,1),data=dexts); #reject normality
mymodel<-garchFit(~arma(0,1) + garch(1,1), cond.dist="QMLE",data=dexts);
summary(mymodel)
plot(mymodel)


#Compare forecast errors

y<-exts
S=round(0.90*length(y))
h=1
error1.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(0,1,1),seasonal=c(0,0,0), method="ML")
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error1.h<-c(error1.h,y[i+h]-predict.h)
}
error2.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(1,0,1),seasonal=c(0,0,0), method="ML")
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error2.h<-c(error2.h,y[i+h]-predict.h)
}
error3.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(1,1,0),seasonal=c(0,0,0), method="ML")
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error3.h<-c(error3.h,y[i+h]-predict.h)
}

cbind(error1.h,error2.h, error3.h)
dm.test(error1.h,error3.h,h=h,power=1) #do not reject, not significantly different from each other
dm.test(error1.h,error2.h,h=h,power=1) #do not reject, not significantly different from each other

MAE1 <- mean(abs(error1.h))
MAE2 <- mean(abs(error2.h))
MAE3 <- mean(abs(error3.h))

MAE1; MAE2; MAE3; 
AIC(arima1, k = log(272))
AIC(arima2, k = log(272))
AIC(arima3, k = log(272))

AIC(arima1)
AIC(arima2)
AIC(arima3)

#Multivariate

#second series

mydata2 <- read.table(file=file.choose(),header=TRUE, dec=".")
attach(mydata2)
names(mydata2)
data <- data.frame(mydata, mydata2)
attach(data)
cpits <- ts(INDCPIALLMINMEI, frequency = 12, start=c(2000,1))
logcpits <- log(cpits)
dcpits<-diff(cpits)
ts.plot(cpits, exts, col=c("black","blue"), ylim=c(30, 150))
summary(cpits)
#plot series in logs and differences

dlogcpits <- diff(diff(logcpits), lag = 12)
plot(cpits)
plot(dcpits)
plot(logcpits)
plot(dlogcpits)
acf(cpits)
acf(dcpits)
acf(logcpits)
acf(dlogcpits)
pacf(cpits)
pacf(dcpits)
pacf(logcpits)
pacf(dlogcpits)
summary(dlogcpits)
monthplot(dlogcpits)
monthplot(dcpits) #seasonality present

#Tests for stationarity and white noise
CADFtest(dcpits, type= "drift", criterion= "BIC", max.lag.y=max.lag) # do not reject, not stationary
Box.test(dlogcpits, lag = max.lag, type = "Ljung-Box") #reject, not white noise
CADFtest(dlogcpits, type= "drift", criterion= "BIC", max.lag.y=max.lag) #reject, stationary
Box.test(dlogcpits, lag = max.lag, type = "Ljung-Box") #reject, not white noise


#cointegration

fit_ci <- lm(exts ~ cpits)
res_fit_ci <- fit_ci$residuals
plot(fit_ci$residuals)
acf(fit_ci$residuals)
CADFtest(res_fit_ci,type="drift",criterion="BIC",max.lag.y=max.lag) #do not reject >0.01, no cointegration
Box.test(fit_ci$residuals, lag = max.lag, type = "Ljung-Box") #reject, not white noise

fit_ci2 <- lm(exts ~ log(cpits))
res_fit_ci2 <- fit_ci2$residuals
plot(fit_ci2$residuals)
acf(fit_ci2$residuals) 
CADFtest(res_fit_ci2,type="drift",criterion="BIC",max.lag.y=max.lag) #do not reject, not stationary 
Box.test(fit_ci2$residuals, lag = max.lag, type = "Ljung-Box") #reject, not white noise



#Granger causality
fit_ci3 <- lm(dexts ~ diff(log(cpits)))
res_fit_ci3 <- fit_ci3$residuals
plot(fit_ci3$residuals)
acf(fit_ci3$residuals) #significant autocorrelation at lag 1
CADFtest(res_fit_ci3,type="drift",criterion="BIC",max.lag.y=max.lag) #reject, stationary
Box.test(fit_ci3$residuals, lag = max.lag, type = "Ljung-Box") #reject, not white noise



#lag=1
fit_adlm1 <- lm(dex.0 ~ dex.1+dlogcpi.1, data=dt)
redfit_adlm1 <- lm(dex.0 ~ dex.1, data=dt)
summary(fit_adlm1)
summary(redfit_adlm1)
acf(fit_adlm1$residuals)
Box.test(fit_adlm1$residuals, lag = max.lag, type = "Ljung-Box") #do not reject, white noise
acf(redfit_adlm1$residuals)
Box.test(redfit_adlm1$residuals, lag = max.lag, type = "Ljung-Box")#do not reject, white noise

anova( fit_adlm1, redfit_adlm1) #do not reject no granger causality, therefore no incremental explanatory power


#VAR model

multivardata <- cbind(dlogcpits, dexts)
multivardata<-na.omit(multivardata)
VARselect(multivardata, lag.max = max.lag, type = "const")
fit_var1<-VAR(multivardata,type="const",p=1)
summary(fit_var1)
var1_residuals<-resid(fit_var1)
par(mfrow=c(2,2))
acf(var1_residuals[,1]) #significant autocorrelation at lag 12
acf(var1_residuals[,2])
ccf(var1_residuals[,1],var1_residuals[,2])
par(mfrow=c(1,1))
Box.test(var1_residuals[,1], lag = max.lag, type = "Ljung-Box") #reject dlogcpi, not white noise
Box.test(var1_residuals[,2], lag = max.lag, type = "Ljung-Box") #do not reject, white noise


fit_var12<-VAR(multivardata,type="const",p=12)
summary(fit_var12)
var12_residuals<-resid(fit_var12)
par(mfrow=c(2,2))
acf(var12_residuals[,1])
acf(var12_residuals[,2])
ccf(var12_residuals[,1],var12_residuals[,2])
par(mfrow=c(1,1))
Box.test(var12_residuals[,1], lag = max.lag, type = "Ljung-Box") #do not reject dlogcpi
Box.test(var12_residuals[,2], lag = max.lag, type = "Ljung-Box") #do not reject
#Impulse-response

irf_var<-irf(fit_var12,ortho=FALSE,boot=TRUE)
plot(irf_var)



#VAR prediction

myforecast <- predict(fit_var12, n.ahead = 12)

par(mfrow=c(2,1))


lastvalue <- rep(c(80.2538), each = 12)
dvarexpected <- predictedvar$fcst$dexts[,1]
dvarexpected<-transform(dvarexpected, cumFreq=cumsum(dvarexpected))
dvarlower <- predictedvar$fcst$dexts[,2]
dvarlower<-transform(dvarlower, cumFreq=cumsum(dvarlower))
dvarupper <- predictedvar$fcst$dexts[,3]
dvarupper<-transform(dvarupper, cumFreq=cumsum(dvarupper))
varexpected <- lastvalue+dvarexpected[,2]
varlower <- lastvalue+dvarlower[,2]
varupper <- lastvalue+dvarupper[,2]
varexpectedts <- ts(varexpected, frequency = 12, start = c(2022,10))
varlowerts <- ts(varlower, frequency = 12, start = c(2022,10))
varupperts <- ts(varupper, frequency = 12, start = c(2022,10))

#Plot the prediction
plot.ts(exts,xlim=c(2020,2025),ylim=c(40,90))
lines(varexpectedts,col="red")
lines(varlowerts,col="blue")
lines(varupperts,col="blue")

#vecm
library(urca)
dt2<-cbind(exts, logcpits)
dt2<-na.omit(dt2)
VARselect(dt2, lag.max = max.lag, type = "const")
trace_test<-ca.jo(dt2,type="trace",K=16,ecdet="const",spec="transitory")
summary(trace_test) #r=1
maxeigen_test<-ca.jo(dt2,type="eigen",K=16,ecdet="const",spec="transitory")
summary(maxeigen_test) #r=1
fit_vecm1<-cajorls(trace_test,r=1)
fit_vecm1
fit_vecm2<-cajorls(maxeigen_test,r=1)
fit_vecm2

fit_var<-vec2var(trace_test,r=1)
myforecast<-predict(fit_var,n.ahead=12)
par(mfrow=c(2,1))

exts_forecast<-ts(myforecast$fcst$exts[,1],frequency=12,start=c(2022,10))
exts_lower<-ts(myforecast$fcst$exts[,2],frequency=12,start=c(2022,10))
exts_upper<-ts(myforecast$fcst$exts[,3],frequency=12,start=c(2022,10))
ts.plot(exts_forecast,exts_lower,exts_upper,col=c("black","red","red"))
title(main = "12-step-ahead forecast of exts")
par(mfrow=c(2,1))


logcpi_forecast<-ts(myforecast$fcst$logcpi[,1],frequency=12,start=c(2022,10))
logcpi_lower<-ts(myforecast$fcst$logcpi[,2],frequency=12,start=c(2022,10))
logcpi_upper<-ts(myforecast$fcst$logcpi[,3],frequency=12,start=c(2022,10))
ts.plot(logcpi_forecast,logcpi_lower,logcpi_upper,col=c("black","red","red"))
title(main = "12-step-ahead forecast of logcpi")

