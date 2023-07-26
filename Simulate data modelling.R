library(dlnm)
library(mgcv)
library(gam)
library(ggplot2)


data("chicagoNMMAPS")
cbCHIpm =crossbasis(chicagoNMMAPS$temp, lag = 5, argvar = list(fun="ns", 3), 
                    arglag = list(fun='ns',6))
CHIglm=gam(death~cbCHIpm:o3,family = quasipoisson(), data = chicagoNMMAPS)
CHIpre=crosspred(cbCHIpm,CHIglm)
plot(CHIpre,ticktype='detailed',xlab="\n\n PM 2.5",ylab="\n\n Lag (weeks)", zlab="\n\n RR",
     border='#3366FF',
     # col='#99FFCC',
     shade = 0.1,
     cex.lab=1.3,cex.axis=1.2,lwd=1.3,
     theta = 40, phi = 28, # theta for left-right, phi for up-down
     ltheta = 35)
