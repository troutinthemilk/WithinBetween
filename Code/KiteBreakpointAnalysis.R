source("Kite_Funcs.R")

#analyze the snail kite example
kite.est	        <- as.matrix(read.csv(file="../Kite/SNKI abundance_south.csv",header=T))
years.kite    		<- kite.est[,1]

Water.new  			<- read.csv(file="../Data/gauges_1996_2016_no_outliers.csv",header=T, sep=',')
WCA3A.dat       <- Water.new$WCA3A
names(WCA3A.dat) <- Water.new[,1]

##fixes NA's in water level data
days    <- 1:length(WCA3A.dat)
lo      <- loess(WCA3A.dat ~ days, span=30/length(WCA3A.dat))
na.vec  <- which(is.na(WCA3A.dat))
lo.predict <- predict(lo, newdata=days, se = FALSE)
WCA3A.dat[na.vec] <- lo.predict[na.vec]

ind         <- 1
water96  		<- WCA3A.dat[-(1:ind)] #series starts on July 1, 1996
water97     <- WCA3A.dat[-(1:(ind+365))] #series starts on July 1, 1997

#get average yearly water level
years.water 	<- levels(years.kite)+1
water.mean		<- vector('numeric',length(years.water))

par(mar=c(5.1,4.1,4.1,4.1))
plot(kite.est, type='b', lwd=2, pch=19, xlab="Year", ylab="Abundance", cex.lab=1.4)
par(new=TRUE)
plot(WCA3A.dat, type="l",col="dodgerblue3",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("Water levels",side=4,line=3, cex=1.4, col="dodgerblue3")
stop()

j <- 1
for(i in 1:(dim(kite.est)[1]+1)) {
  water.mean[j] <- mean(WCA3A.dat[((i-1)*365+1):(365*i)])
	j <- j + 1
}

kite.lag    <- log(kite.est[1:(length(years.kite)-1),2])
kite.pgr    <- log(kite.est[2:(length(years.kite)),2])#/(kite.est[1:(length(years.kite)-1),2]))

#running calculation of model fits
if(F) {
  ll.vec <- vector('numeric', 16)
  cor.vec <- ll.vec
  for(i in 5:16) {
    mean.fits		    <- lm.fits(kite.pgr=kite.pgr[1:i], kite.lag=kite.lag[1:i], water.cov=water.mean[1:i+1])
    ll.vec[i] <- mean.fits[[1]][2]/i
    cor.vec[i] <- cor(kite.pgr[1:i], water.mean[2:(i+1)])
  }

  plot(years.kite[-(1:5)], exp(ll.vec[5:16]), type='l', lwd=2, xlab="Cumulative years considered", ylab="Average likelihood of an observation")
  plot(years.kite[-(1:5)], cor.vec[5:16], type='l', lwd=2, xlab="Cumulative years considered", ylab="Correlation between pgr and mean water level")
}

#fit model with breakpoint in tau and d.
lm.mod    <- lm(kite.pgr ~ kite.lag)
lm.mod1   <- lm(kite.pgr[1:10] ~ kite.lag[1:10])
lm.mod2   <- lm(kite.pgr[11:16] ~ kite.lag[11:16])

par.vec             <- c(a=coef(lm.mod)[1], b=log(abs(coef(lm.mod)[2])), sigma=log(sd(kite.pgr)))
breakpoint_nocov    <- mlsl(x0=par.vec, fn=ricker.generic, lower=c(a=-10, b=-10, sig=-5), upper=c(a=10, b=2, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10)
names(breakpoint_nocov$par) <- names(par.vec)


par.vec   <- c(a1=coef(lm.mod)[1], a2=coef(lm.mod)[1], b=log(abs(coef(lm.mod)[2])), sigma=log(sd(kite.pgr)))
breakpoint_nocov_a   <- mlsl(x0=par.vec, fn=ricker.generic, lower=c(a1=-10, a2=-10, b=-10, sig=-5), upper=c(a1=10, a2=10, b=2, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10, a.flag=T)
names(breakpoint_nocov_a$par) <- names(par.vec)

par.vec   <- c(a=coef(lm.mod)[1], b1=log(abs(coef(lm.mod)[2])), b2=log(abs(coef(lm.mod)[2])), sigma=log(sd(kite.pgr)))
breakpoint_nocov_b   <- mlsl(x0=par.vec, fn=ricker.generic, lower=c(a=-10, b1=-10, b1=-10, sig=-5), upper=c(a=10, b1=2, b2=2, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10, b.flag=T)
names(breakpoint_nocov_b$par) <- names(par.vec)

par.vec   <- c(a1=coef(lm.mod)[1], a2=coef(lm.mod)[1], b1=log(abs(coef(lm.mod)[2])), b2=log(abs(coef(lm.mod)[2])), sigma=log(sd(kite.pgr)))
breakpoint_nocov_ab   <- mlsl(x0=par.vec, fn=ricker.generic, lower=c(a1=-20, a2=-20, b1=-100, b1=-100, sig=-5), upper=c(a1=50, a2=50, b1=10, b2=10, sig=10), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e4), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10, a.flag=T, b.flag=T)
names(breakpoint_nocov_ab$par) <- names(par.vec)

#stop()
ll.vec <- c(breakpoint_nocov$value, breakpoint_nocov_a$value, breakpoint_nocov_b$value, breakpoint_nocov_ab$value)
#11.382879  7.916329  7.205119  4.140605

n <- length(kite.pgr)
k <- c(length(breakpoint_nocov$par), length(breakpoint_nocov_a$par), length(breakpoint_nocov_b$par), length(breakpoint_nocov_ab$par))
print(BIC.cus(ll.vec, k, n))

#print(AICc.cus(ll.vec, k, n))
# 30.76576 27.46902 26.04660 24.28121

library(numDeriv)
hmat <- hessian(func=ricker.generic, x=breakpoint_nocov_ab$par, kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10, a.flag=T, b.flag=T, method.args=list(eps=1e-6, d=0.001))
breakpointse <- solve(hmat)

stop()

par.vec   <- c(a1=coef(lm.mod)[1], a2=coef(lm.mod)[1], b1=log(abs(coef(lm.mod)[2])), b2=log(abs(coef(lm.mod)[2])), d1=0, sigma=log(sd(kite.pgr)))
breakpoint_disc_ab1   <- mlsl(x0=par.vec, fn=ricker.generic.discov1, lower=c(a1=-10, a2=-10, b1=-10, b1=-10, d1=-10, sig=-5), upper=c(a1=10, a2=10, b1=2, b2=2, d1=10, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=water.mean, bp=10, a.flag=T, b.flag=T)
names(breakpoint_disc_ab1$par) <- names(par.vec)

par.vec   <- c(a1=coef(lm.mod)[1], a2=coef(lm.mod)[1], b1=log(abs(coef(lm.mod)[2])), b2=log(abs(coef(lm.mod)[2])), d1=0, sigma=log(sd(kite.pgr)))
breakpoint_disc_ab2   <- mlsl(x0=par.vec, fn=ricker.generic.discov2, lower=c(a1=-10, a2=-10, b1=-10, b1=-10, d1=-10, sig=-5), upper=c(a1=10, a2=10, b1=2, b2=2, d1=10, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=water.mean[-1], bp=10, a.flag=T, b.flag=T)
names(breakpoint_disc_ab2$par) <- names(par.vec)


ll.vec <- c(breakpoint_disc_ab1$value, breakpoint_disc_ab2$value)
n <- length(kite.pgr)
k <- c(length(breakpoint_disc_ab1$par), length(breakpoint_disc_ab2$par))
print(BIC.cus(ll.vec, k, n))
#print(AICc.cus(ll.vec, k, n))

par.vec   <- c(a1=coef(lm.mod)[1], a2=coef(lm.mod)[1], b1=log(abs(coef(lm.mod)[2])), b2=log(abs(coef(lm.mod)[2])), d1=0, delta=log(180), sigma=log(sd(kite.pgr)))
breakpoint_cont_ab1   <- mlsl(x0=par.vec, fn=ricker.contcov1, lower=c(a1=-10, a2=-10, b1=-10, b1=-10, d1=-10, delta=-log(1), sig=-5), upper=c(a1=10, a2=10, b1=2, b2=2, d1=10, delta=log(5*365), sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat[-(1:(181))], tau1=0, tau2=0, bp=10, a.flag=T, b.flag=T)
names(breakpoint_cont_ab1$par) <- names(par.vec)

par.vec   <- c(a1=coef(lm.mod)[1], a2=coef(lm.mod)[1], b1=log(abs(coef(lm.mod)[2])), b2=log(abs(coef(lm.mod)[2])), d1=0, delta=log(180), sigma=log(sd(kite.pgr)))
breakpoint_cont_ab2   <- mlsl(x0=par.vec, fn=ricker.contcov2, lower=c(a1=-10, a2=-10, b1=-10, b1=-10, d1=-10, delta=log(1), sig=-5), upper=c(a1=10, a2=10, b1=2, b2=2, d1=10, delta=log(2*365), sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, tau1=182, tau2=182, bp=10, a.flag=T, b.flag=T)
names(breakpoint_cont_ab2$par) <- names(par.vec)

#breakpoint_cont_ab2   <- optim(par=breakpoint_cont_ab2$par, fn=ricker.contcov2, local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, method="BFGS", kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, tau1=182, tau2=182, bp=10, a.flag=T, b.flag=T)

library(numDeriv)
hmat <- hessian(func=ricker.contcov2, x=breakpoint_cont_ab2$par, kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, tau1=182, tau2=182, bp=10, a.flag=T, b.flag=T, method.args=list(eps=1e-4, d=0.00001))
breakpointse <- solve(hmat)

ll.vec <- c(breakpoint_cont_ab1$value, breakpoint_cont_ab2$value)
n <- length(kite.pgr)
k <- c(length(breakpoint_cont_ab1$par), length(breakpoint_cont_ab2$par))
print(BIC.cus(ll.vec, k, n))
#print(AICc.cus(ll.vec, k, n))

stop()
par.vec   <- c(a=coef(lm.mod)[1], b=log(abs(coef(lm.mod)[2])), d=0, delta=log(180), tau=10, sigma=log(sd(kite.pgr)))
breakpoint_null   <- mlsl(x0=par.vec, fn=ricker.decay.generic, lower=c(a=-10, b=-10, d=0, delta=log(1), tau=1, sig=-5), upper=c(a=10, b=2, d=10, delta=log(2*365), tau=366+181, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10)
names(breakpoint_null$par) <- names(par.vec)

par.vec   <- c(a=coef(lm.mod)[1], b1=log(abs(coef(lm.mod)[2])), b2=log(abs(coef(lm.mod)[2])), d=0, delta=log(180), tau=10, sigma=log(sd(kite.pgr)))
breakpoint_b   <- mlsl(x0=par.vec, fn=ricker.decay.generic, lower=c(a=-10, b1=-10, b2=-10, d=0, delta=log(1), tau=1, sig=-5), upper=c(a=10, b1=2, b2=2, d=10, delta=log(2*365), tau=366+181, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10, b.flag=T)
names(breakpoint_b$par) <- names(par.vec)
#wiht 2 b's: 6.437419
#with 1 b:  10.65091

par.vec   <- c(a1=coef(lm.mod)[1], a2=coef(lm.mod)[1], b=log(abs(coef(lm.mod)[2])), d=0, delta=log(180), tau=10, sigma=log(sd(kite.pgr)))
breakpoint_a   <- mlsl(x0=par.vec, fn=ricker.decay.generic, lower=c(a1=-10, a2=-10, b=-10, d=0, delta=log(1), tau=1, sig=-5), upper=c(a1=10, a2=10, b=2, d=10, delta=log(2*365), tau=366+181, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10, a.flag=T)
names(breakpoint_a$par) <- names(par.vec)

par.vec   <- c(a1=coef(lm.mod)[1], a2=coef(lm.mod)[1], b1=log(abs(coef(lm.mod)[2])), b2=log(abs(coef(lm.mod)[2])), d=0, delta=log(180), tau=10, sigma=log(sd(kite.pgr)))
breakpoint_ab   <- mlsl(x0=par.vec, fn=ricker.decay.generic, lower=c(a1=-10, a2=-10, b1=-10, b2=-10, d=0, delta=log(1), tau=1, sig=-5), upper=c(a1=10, a2=10, b1=2, b2=2, d=10, delta=log(2*365), tau=366+181, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10, a.flag=T, b.flag=T)
names(breakpoint_ab$par) <- names(par.vec)


par.vec   <- c(a=coef(lm.mod)[1], b=log(abs(coef(lm.mod)[2])), d=0, delta=log(180), tau1=10, tau2=10, sigma=log(sd(kite.pgr)))
breakpoint_tau   <- mlsl(x0=par.vec, fn=ricker.decay.generic, lower=c(a=-10, b=-10, d1=0, delta=log(1), tau1=0, tau2=0, sig=-5), upper=c(a=10, b1=2, d1=10, delta=log(2*365), tau1=366+181, tau1=366+181, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10, tau.flag=T)
names(breakpoint_tau$par) <- names(par.vec)

par.vec   <- c(a=coef(lm.mod)[1], b=log(abs(coef(lm.mod)[2])), d1=0, d2=0, delta=log(180), tau1=10, tau2=10, sigma=log(sd(kite.pgr)))
breakpoint_tau_delta   <- mlsl(x0=par.vec, fn=ricker.decay.generic, lower=c(a=-10, b=-10, d1=0, delta1=log(1), delta2=log(1), tau1=0, tau2=0, sig=-5), upper=c(a=10, b1=2, d=10, delta1=log(2*365), delta2=log(2*365), tau1=366+181, tau2=366+181, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10, tau.flag=T, delta.flag=T)
names(breakpoint_tau_delta$par) <- names(par.vec)

par.vec   <- c(a=coef(lm.mod)[1], b=log(abs(coef(lm.mod)[2])), d=0, delta1=log(180), delta2=log(180), tau1=10, sigma=log(sd(kite.pgr)))
breakpoint_delta   <- mlsl(x0=par.vec, fn=ricker.decay.generic, lower=c(a=-10, b=-10, d=0, delta1=log(1), delta2=log(1), tau=1, sig=-5), upper=c(a=10, b1=2, d=10, delta1=log(2*365), delta2=log(2*365), tau=366+181, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10, delta.flag=T)
names(breakpoint_delta$par) <- names(par.vec)

par.vec   <- c(a=coef(lm.mod)[1], b=log(abs(coef(lm.mod)[2])), d1=0, d2=0, delta1=log(180), delta2=log(180), tau1=10, sigma=log(sd(kite.pgr)))
breakpoint_delta_d   <- mlsl(x0=par.vec, fn=ricker.decay.generic, lower=c(a=-10, b=-10, d1=0, d2=0, delta1=log(1), delta2=log(1), tau=1, sig=-5), upper=c(a=10, b1=2, d1=10, d2=10, delta1=log(2*365), delta2=log(2*365), tau=366+181, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=WCA3A.dat, bp=10, d.flag=TRUE, delta.flag=TRUE)
names(breakpoint_delta_d$par) <- names(par.vec)

#breakpoint_tau$value=4.152476
ll.vec <- c(breakpoint_null$value, breakpoint_a$value, breakpoint_ab$value, breakpoint_b$value, breakpoint_tau$value, breakpoint_tau_delta$value, breakpoint_delta$value, breakpoint_delta_d$value)
#9.076348 5.628098 6.437419 8.165949 6.383364 7.133756 6.212037

n <- length(kite.pgr)
k <- c(length(breakpoint_null$par), length(breakpoint_a$par), length(breakpoint_ab$par), length(breakpoint_b$par), length(breakpoint_tau$par), length(breakpoint_tau_delta$par), length(breakpoint_delta$par), length(breakpoint_delta_d$par))
print(BIC.cus(ll.vec, k, n))
#print(AICc.cus(ll.vec, k, n))

save.image(file="KiteBreakpointModels.Rdata")

