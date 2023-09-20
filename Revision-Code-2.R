
library(mgcv)
library(MASS)
library(dlnm)
library(DHARMa)
library(car)
library(ggplot2)
library(splines) 
library(tsModel)

attach(sh_sum)

qaic <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum( dpois( model$y, model$fitted.values, log=TRUE) )
  return(-2*loglik + 2*dim(as.matrix(model[["coefficients"]]))[1]*phi)
}

# Univariable -------------------------------------------------------------


# PM25 ###############################################################################

plot(pm25, rate)
plot(smooth.spline(pm25,rate))
varknots=c(50,100,150)
varknots=equalknots(pm25,fun = 'ns',3)
basisPM=crossbasis(sh_sum$pm25, lag = 5, 
                   argvar = list(fun="ns", 4, knots=c(50,150)),
                   arglag = list(fun='ns',6, knots=c(2))) # cb.pm25

## df on var and lag do not change qaic in the model 

model0=      gam(posi~ basisPM+ns(time,8)+ns(week,6)+offset(log(case)), 
                 family = quasipoisson(link = 'log'), sh_sum)

qaic(model0)
summary(model0)

pred1 <- crosspred(basisPM,model0,cen = median(sh_sum$pm25))

rr.ci.lag(pred1,3) # first max
# 1.319 (1.041 1.672) var=155
rr.ci.lag(pred1,6) # highest
# 1.467 (0.762 2.824) var=230
rr.ci.cum(pred1) # overall max
# 3.404 (1.532 7.560) var=165



par(mgp=c(2.5,0.5,0))
plot(pred1, "contour", xlab="pm25", key.title=title("RR"),cex.axis=2,
     plot.axes={axis(1,cex.axis=1)
       axis(2,cex.axis=1)},
     key.axes = axis(4,cex.axis=1),
     plot.title=title(xlab=" PM 2.5", ylab="Lag (weeks)",cex.main=2,cex.lab=1.5))

plot(pred1,
     ticktype='detailed', #刻度类型
     ntick=4, #刻度数量
     border='#3366FF',
     xlab="\n\n PM 2.5",ylab="\n\n Lag (weeks)", zlab="\n\n RR",
     # col='#99FFCC',
     shade = 0.1,
     cex.lab=1.3,cex.axis=01,lwd=1.3,
     theta = 40, phi = 28, # theta for left-right, phi for up-down
     ltheta = 35)

plot(pred1,'slices',var=180,lwd=2,ci.arg=list(density=20,col=2),
     main="Lag-response for different AH levels", col=2, xlab='PM 2.5', ylab='RR')
plot(pred1,'slices',lag=1,lwd=2,ci.arg=list(lty=3,col='grey'),
     main="Lag-response for different AH levels", col=2, xlab='PM 2.5', ylab='RR')

plot(pred1,'overall',lwd=2,ci.arg=list(density=20,col='grey'),
     main="Lag-response for different AH levels", col=2, xlab='PM 2.5', ylab='RR')


plot(pred1,'overall', ci.arg = list(lty=3,col='grey'),

     xlab=expression(paste(PM[2.5], " (", mu, g, "/", m^3, ")", sep = "")), ylab='RR', 
     lwd=3,col='red')

plot(pre.pm25,'slices',lag=1, ci.arg = list(lty=3,col='grey'), xlab=expression(paste(PM[2.5], " (", mu, g, "/", m^3, ")", sep = "")), ylab='RR', lwd=3,col='red')
plot(pred1,'slices',lag=1, ci.arg = list(lty=3,col='grey'), xlab=expression(paste(PM[2.5], " (", mu, g, "/", m^3, ")", sep = "")), ylab='RR', lwd=3,col='red')


# AH  ##########################################################################



plot(smooth.spline(ah,rate))

varknots <- equalknots(ah,fun="ns",df=6)

basisAH=crossbasis(sh_sum$ah, lag = 5, 
                   argvar = list(fun="ns", 4, knots=c(10,15)),
                   arglag = list(fun='ns',6, knots=c(2))) # cb.ah


model0=      gam(posi~ basisAH +ns(time,8)+ns(week,6)+offset(log(case)), 
                 family = quasipoisson(link = 'log'), sh_sum)

qaic(model0)

pred1 <- crosspred(basisAH,model0,cen = median(sh_sum$ah))

rr.ci.lag(pred1,4) # highest
# 2.001 (1.278 3.162) var=4

rr.ci.cum(pred1) # overall max
# 31.534 (5.021 198.043) var=4

par(mgp=c(2.5,0.5,0))
plot(pred1, "contour", xlab="AH", key.title=title("RR"),cex.axis=2,
     plot.axes={axis(1,cex.axis=1)
       axis(2,cex.axis=1)},
     key.axes = axis(4,cex.axis=1),
     plot.title=title(xlab=" AH", ylab="Lag (weeks)",cex.main=2,cex.lab=1.5))

plot(pred1,
     ticktype='detailed', #刻度类型
     ntick=4, #刻度数量
     border='#3366FF',
     xlab="\n\n AH",ylab="\n\n Lag (weeks)", zlab="\n\n RR",
     # col='#99FFCC',
     shade = 0.1,
     cex.lab=1.3,cex.axis=01,lwd=1.3,
     theta = 40, phi = 28, # theta for left-right, phi for up-down
     ltheta = 35)

plot(pred1,'overall',lwd=2,ci.arg=list(density=20,col=2), ylim=c(0,30),
     main="Lag-response for different AH levels", col=2, xlab='AH', ylab='RR')

plot(pred1,'overall', ci.arg = list(lty=3,col='grey'),
     ylim=c(0,30),
     xlab='AH', ylab='RR', 
     lwd=3,col='red')
# Min Temp #####################################################################

plot(tempmin, rate)
smooth.spline(tempmin,rate)
plot(smooth.spline(tempmin,rate))

basisTEMP=crossbasis(sh_sum$tempmin, lag = 5, 
                   argvar = list(fun="ns", 4, knots=c(5,15,20)),
                   arglag = list(fun='ns',6, knots=c(2))) # cb.tempmin


model0=      gam(posi~ basisTEMP +ns(time,8)+ns(week,6)+offset(log(case)), 
                 family = quasipoisson(link = 'log'), sh_sum)

qaic(model0)
# summary(model0)
pred1 <- crosspred(basisTEMP,model0,cen = median(sh_sum$tempmin))


rr.ci.lag(pred1,3) # first max
# 2.188 (1.414 3.387) var=0.5
rr.ci.lag(pred1,6) # highest
# 1.251 (0.624 2.505) var=22
rr.ci.lag(pred1,1)
# 1.271 (1.032 1.565)
rr.ci.cum(pred1) # overall max
# 18.836 (4.366 81.271) var=2

par(mgp=c(2.5,0.5,0))
plot(pred1, "contour", xlab="Min Temp", key.title=title("RR"),cex.axis=2,
     plot.axes={axis(1,cex.axis=1)
       axis(2,cex.axis=1)},
     key.axes = axis(4,cex.axis=1),
     plot.title=title(xlab=" Min Temp", ylab="Lag (weeks)",cex.main=2,cex.lab=1.5))

plot(pred1,
     ticktype='detailed', #刻度类型
     ntick=4, #刻度数量
     border='#3366FF',
     xlab="\n\n Min Temp",ylab="\n\n Lag (weeks)", zlab="\n\n RR",
     # col='#99FFCC',
     shade = 0.1,
     cex.lab=1.3,cex.axis=01,lwd=1.3,
     theta = 40, phi = 28, # theta for left-right, phi for up-down
     ltheta = 35)

plot(pred1,'overall',lwd=2,ci.arg=list(density=20,col=2),
     main="Lag-response for different AH levels", col=2, xlab='Min Temp', ylab='RR')

plot(pred1,'overall', ci.arg = list(lty=3,col='grey'),
     ylim=c(0,30),
     xlab='Min Temp', ylab='RR', 
     lwd=3,col='red')




# RH ####
smooth.spline(tempmin,rate)
plot(smooth.spline(rh,rate))

basisRH=crossbasis(sh_sum$rh, lag = 5, 
                     argvar = list(fun="ns", 4),
                     arglag = list(fun='ns',6, knots=c(2))) # cb.tempmin


model0=      gam(posi~ basisRH +ns(time,8)+ns(week,6)+offset(log(case)), 
                 family = quasipoisson(link = 'log'), sh_sum)

qaic(model0)
# summary(model0)
pred1 <- crosspred(basisRH,model0,cen = median(sh_sum$tempmin))

plot(pred1, "contour", xlab="Min Temp", key.title=title("RR"),cex.axis=2,
     plot.axes={axis(1,cex.axis=1)
       axis(2,cex.axis=1)},
     key.axes = axis(4,cex.axis=1),
     plot.title=title(xlab=" Min Temp", ylab="Lag (weeks)",cex.main=2,cex.lab=1.5))

plot(pred1,
     ticktype='detailed', #刻度类型
     ntick=4, #刻度数量
     border='#3366FF',
     xlab="\n\n RH",ylab="\n\n Lag (weeks)", zlab="\n\n RR",
     # col='#99FFCC',
     shade = 0.1,
     cex.lab=1.3,cex.axis=01,lwd=1.3,
     theta = 40, phi = 28, # theta for left-right, phi for up-down
     ltheta = 35)


# Interactive effect -----------------------------------------------------



#### AH Interaction (Cen: 5,15)  ####

basisPM=crossbasis(sh_sum$pm25, lag = 5, 
                   argvar = list(fun="ns", 4, knots=c(50,150)),
                   arglag = list(fun='ns',6,knots=c(2))) # cb.pm25

model.AH0=      gam(posi~ basisPM +ah+offset(log(case))+ns(time,8)+ns(week,6), 
                 family = quasipoisson(link = 'log'), sh_sum)

qaic(model.AH0)

AH.mov.avg=runMean(ah)
AHc3=AH.mov.avg -5
AHc8=AH.mov.avg -15

basis1=basisPM*AHc3
basis2=basisPM*AHc8

model.AH1 = update(model.AH0, .~. + basis1)
model.AH2 = update(model.AH0, .~. + basis2)

pred1.ah.int =  crosspred(basisPM,model.AH1,bylag = 0.5,cen=180)
pred2.ah.int =  crosspred(basisPM,model.AH2,bylag = 0.5,cen=180)

rr.ci.lag(pred1.ah.int,1)
# $rr
# [1] 5.633583
# 
# $ci.low
# [1] 0.7689219
# 
# $ci.high
# [1] 41.27501
# 
# $lag
# [1] 0
# 
# $var
# [1] "230"
rr.ci.lag(pred1.ah.int,6)
# 1.944 (0.191 19.789)
rr.ci.lag(pred2.ah.int,6)
# 1.897 (0.617 5.830)
## Plots

par(mgp=c(2.5,0.5,0))

plot(pred1.ah.int,
     ticktype='detailed', #刻度类型
     ntick=4, #刻度数量
     border='#3366FF',
     xlab="\n\n PM 2.5 and AH",ylab="\n\n Lag (weeks)", zlab="\n\n RR",
     # col='#99FFCC',
     shade = 0.1,
     cex.lab=1.3,cex.axis=01,lwd=1.3,
     theta = 40, phi = 28, # theta for left-right, phi for up-down
     ltheta = 35)

plot(pred2.ah.int,
     ticktype='detailed', #刻度类型
     ntick=4, #刻度数量
     border='#3366FF',
     xlab="\n\n PM 2.5 and AH",ylab="\n\n Lag (weeks)", zlab="\n\n RR",
     # col='#99FFCC',
     shade = 0.1,
     cex.lab=1.3,cex.axis=01,lwd=1.3,
     theta = 40, phi = 28, # theta for left-right, phi for up-down
     ltheta = 35)


par(mfrow=c(1,2))

plot(pred1.ah.int,var=80,lwd=2,ci.arg=list(density=20,col=2),ylim=c(0,5),
     #main="Lag-response for different AH", 
     col=2, ylab='RR')
lines(pred2.ah.int,var=80,ci="area",lwd=2,ci.arg=list(density=20,col=4,angle=-45),col=4)
legend("top",c('PM 2.5=80',"AH = 5","AH = 15"),lwd=2,
       col=c(0,2,4),inset=0.2)

plot(pred1.ah.int,var=150,lwd=2,ci.arg=list(density=20,col=2),ylim=c(-1,5),
     main="Lag-response for different AH levels", col=2, ylab='RR')
lines(pred2.ah.int,var=150,ci="area",lwd=2,ci.arg=list(density=20,col=4,angle=-45),col=4)
legend("top",c("AH equal to 5","AH equal to 15"),lwd=2,
       col=c(2,4),inset=0.2)

plot(pred1.ah.int,var=c(200),lwd=2,ci.arg=list(density=20,col=2),ylim=c(0,5),
     #main="Lag-response for different AH", 
     col=2, ylab='RR')
lines(pred2.ah.int,var=c(200),ci="area",lwd=2,ci.arg=list(density=20,col=4,angle=-45),col=4)
legend("top",c('PM 2.5=200',"AH = 5","AH = 15"),lwd=2,
       col=c(0,2,4),inset=0.2)

plot(pred1.ah.int,'overall',lwd=2,ci.arg=list(density=20,col=2),ylim=c(0,5),
     main="Lag-response for different AH levels", col=2, ylab='RR')
lines(pred2.ah.int,lag = 4,ci="area",lwd=2,ci.arg=list(density=20,col=4,angle=-45),col=4)
legend("top",c("AH = 5","AH = 15"),lwd=2,
       col=c(2,4),inset=0.2)


par(mfow=c(2,2))
plot(pred1,'slices',var=200, 
     ci.arg = list(lty=3,col='grey'), 
     xlab='PM 2.5 & AH', ylab='RR', lwd=3,col='red')




#### Min Temp Interaction (Cen: 5,20) ####

basisPM=crossbasis(sh_sum$pm25, lag = 5, 
                   argvar = list(fun="ns", 4, knots=c(50,150)),
                   arglag = list(fun='ns',6, knots=c(2))) # cb.pm25

model0=      gam(posi~ basisPM +offset(log(case))+ns(time,8)+ns(week,6), 
                 family = quasipoisson(link = 'log'), sh_sum)
qaic(model0)

summary(sh_sum$tempmin)


AHc3=sh_sum$tempmin -5
AHc8=sh_sum$tempmin -20



baisi0=basisPM
basis1=basisPM*AHc3
basis2=basisPM*AHc8

model_origin = update(model0, .~. + basis0)
model1 = update(model0, .~. + basis1)
model2 = update(model0, .~. + basis2)


pred1 <- crosspred(basisPM,model1,cen = 180)
pred2 <- crosspred(basisPM,model2,cen = 180)

rr.ci.lag(pred1,6)
# 3.129 (0.185 53.001)
rr.ci.lag(pred2,6)
# 4.156 (0.943 18.303)
## Plots

rr.ci(max,pred1)

par(mgp=c(2.5,0.5,0))

plot(pred1,
     ticktype='detailed', #刻度类型
     ntick=4, #刻度数量
     border='#3366FF',
     xlab="\n\n PM 2.5 and Min Temp",ylab="\n\n Lag (weeks)", zlab="\n\n RR",
     # col='#99FFCC',
     shade = 0.1,
    
     cex.lab=1.3,cex.axis=01,lwd=1.3,
     theta = 40, phi = 28, # theta for left-right, phi for up-down
     ltheta = 35)

plot(pred2,
     ticktype='detailed', #刻度类型
     ntick=4, #刻度数量
     border='#3366FF',
     xlab="\n\n PM 2.5 and Min Temp",ylab="\n\n Lag (weeks)", zlab="\n\n RR",
     # col='#99FFCC',
     shade = 0.1,

     cex.lab=1.3,cex.axis=01,lwd=1.3,
     theta = 40, phi = 28, # theta for left-right, phi for up-down
     ltheta = 35)

plot(pred1,var=c(80),lwd=2,ci.arg=list(density=20,col=2),ylim=c(0,5),
     #main="Lag-response for different AH", 
     col=2, ylab='RR')
lines(pred2,var=c(80),ci="area",lwd=2,ci.arg=list(density=20,col=4,angle=-45),col=4)
legend("top",c('PM 2.5=80',"AH = Min Temp=5","Min Temp=20"),lwd=2,
       col=c(0,2,4),inset=0.2)

plot(pred1,var=c(200),lwd=2,ci.arg=list(density=20,col=2),ylim=c(0,5),
     #main="Lag-response for different AH", 
     col=2, ylab='RR')
lines(pred2,var=c(200),ci="area",lwd=2,ci.arg=list(density=20,col=4,angle=-45),col=4)
legend("top",c('PM 2.5=200',"Min Temp=5","Min Temp=20"),lwd=2,
       col=c(0,2,4),inset=0.2)


# RH interactive #### 
basisPM=crossbasis(sh_sum$pm25, lag = 5, 
                   argvar = list(fun="ns", 4, knots=c(50,150)),
                   arglag = list(fun='ns',6, knots=c(2))) # cb.pm25

model0=      gam(posi~ basisPM +offset(log(case))+ns(time,8)+ns(week,6), 
                 family = quasipoisson(link = 'log'), sh_sum)
qaic(model0)

summary(sh_sum$tempmin)


AHc3=sh_sum$tempmin -5
AHc8=sh_sum$tempmin -20



baisi0=basisPM
basis1=basisPM*AHc3
basis2=basisPM*AHc8

model_origin = update(model0, .~. + basis0)
model1 = update(model0, .~. + basis1)
model2 = update(model0, .~. + basis2)


pred1 <- crosspred(basisPM,model1,cen = 180)
pred2 <- crosspred(basisPM,model2,cen = 180)

rr.ci.lag(pred1,6)
# 3.129 (0.185 53.001)
rr.ci.lag(pred2,6)
# 4.156 (0.943 18.303)
## Plots