library(dlnm)
library(gam)
library(car)
library(mgcv)
library(ggplot2)
library(DHARMa)

### The dataset in this study named sh_sum

# Univariable GAMs --------------------------------------------------------
model=gam(posi~pm25*ah+ns(week,6)+ns(time,8)+offset(log(cases)),
          family = poisson(link = 'log'), sh_sum)
summary(model)
simulateResiduals(model)
AIC(model)
sum(residuals(model, type = "pearson")^2) / df.residual(model) # dispersion parameter

# Univariable DLNM --------------------------------------------------------
cb.pm25 =crossbasis(sh_sum$pm25, lag = 5, 
                    argvar = list(fun="ns", 3, knots=c(50,100,150)), arglag = list(fun='ns',6, knots=c(2)))

model=gam(posi ~ cb.pm25+ns(time,8)+ns(week,6)+offset(log(case)), 
          family = poisson(link = 'log'), df)

pre.pm25=crosspred(cb.pm25, model, cen = round(median(sh_sum$rh)),cumul = T)

# 2d heat plot
par(mgp=c(2.3,1,0))
plot(pre.pm25, "contour", xlab="PM 2.5",
     key.title=title('RR'),cex.axis=0.8,
     
     plot.axes={axis(1,cex.axis=1.4)
       axis(2,cex.axis=1.4)},
     
     key.axes = axis(4,cex.axis=1.5), # right side bar plot, cex for text size 
     plot.title=title(xlab="PM 2.5",ylab="Lag (weeks)",cex.main=2,cex.lab=1.5))
# 3d surface approach plot
plot(pre.pm25,ticktype='detailed',xlab="\n\n PM 2.5",ylab="\n\n Lag (weeks)", zlab="\n\n RR",
     border='#3366FF',
     # col='#99FFCC',
     shade = 0.1,
     cex.lab=1.3,cex.axis=1.2,lwd=1.3,
     theta = 40, phi = 28, # theta for left-right, phi for up-down
     ltheta = 35)
# expo-res curve overall RR
plot(pre.pm25,'overall', ci.arg = list(lty=3,col='grey'), xlab=expression(paste(PM[2.5], " (", mu, g, "/", m^3, ")", sep = "")), ylab='RR', lwd=3,col='red')



# Interactive effect DLNM -------------------------------------------------

# PM 2.5 & AH 
int.model.ah=gam(posi~cb.pm25:ah +ns(time,8)+ns(week,6)+offset(log(case)), 
                 family = poisson(link = 'log'), sh_sum)   

pre.pm25.int.ah=crosspred(cb.pm25, int.model.ah,
                          cen = median(sh_sum$ah)
)
par(mgp=c(2.5,0.5,0))
# 2D heat plot
plot(pre.pm25.int.ah, "contour", xlab="pm25", key.title=title("RR"),cex.axis=2,
     plot.axes={axis(1,cex.axis=1)
       axis(2,cex.axis=1)},
     key.axes = axis(4,cex.axis=1),
     plot.title=title(xlab=" PM 2.5 and AH", ylab="Lag (weeks)",cex.main=2,cex.lab=1.5))

# 3D surface plot
plot(pre.pm25.int.ah,
     ticktype='detailed', #刻度类型
     ntick=4, #刻度数量
     border='#3366FF',
     xlab="\n\n PM 2.5 and AH",ylab="\n\n Lag (weeks)", zlab="\n\n RR",
     # col='#99FFCC',
     shade = 0.1,
     cex.lab=1.3,cex.axis=01,lwd=1.3,
     theta = 40, phi = 28, # theta for left-right, phi for up-down
     ltheta = 35)

plot(pre.pm25.int.ah,'overall', ci.arg = list(lty=3,col='grey'), xlab='PM 2.5 & AH', ylab='RR', lwd=3,col='red')


