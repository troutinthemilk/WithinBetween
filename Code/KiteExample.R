library(AICcmodavg)

#fits of all variations on the linear model
lm.fits	<- function(kite.pgr, kite.lag, water.cov, water.rmvec=NULL) {

  ind <- 1:length(kite.lag)
	if(!is.null(water.rmvec)) {water.cov <- water.cov[water.rmvec]}

	water.cov 	<- water.cov[ind]#(water.cov[ind]-mean(water.cov[ind]))/sd(water.cov[ind])

	ricker0.lm  <- lm(kite.pgr ~ 1)
	ricker1.lm  <- lm(kite.pgr ~ kite.lag)
	ricker2.lm  <- lm(kite.pgr ~ kite.lag + water.cov)

#print(names(summary(ricker0.lm)))
#stop()
	aicc.table <- aictab(list(ricker0.lm, ricker1.lm, ricker2.lm), modnames=c('Intercept','DD','DD+Cov'), second.ord=T)
  rsq.table <- c(summary(ricker0.lm)$r.squared, summary(ricker1.lm)$r.squared, summary(ricker2.lm)$r.squared)

  LL.table <- c(logLik(ricker0.lm), logLik(ricker1.lm), logLik(ricker2.lm))
  names(LL.table) <- c('Intercept','DD','DD+Cov')
    #print(LL.table)
    print(summary(ricker2.lm))
	BIC.table <- c(BIC(ricker0.lm), BIC(ricker1.lm), BIC(ricker2.lm))
	names(BIC.table) <- c('Intercept','DD','DD+Cov')
	
  return(BIC.table)
}


#fit the ricker model to discretely sampled data
mismatch.optim <- function(kite.pgr, kite.lag, water.cov) {

  	delta.vec <- 1/seq(7, 365*2, by=1)
  	llVec 	<- vector('numeric', length(delta.vec))
  	llList 	<- list()
  	ll2Vec 	<- vector('numeric', length(delta.vec))
  	ll2List <- list()
  
  	for(i in 1:length(delta.vec)) {
    	llList[[i]] <- ricker.decay(delta.val=delta.vec[i], kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=water.cov)
    	llVec[i]    <- logLik(llList[[i]])[1]

    	ll2List[[i]] <- ricker.decaysq(delta.val=delta.vec[i], kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=water.cov)
    	ll2Vec[i]    <- logLik(llList[[i]])[1]

  	}

  mBest 	 	<- llList[[which.max(llVec)]]
  k 		 	<- length(coef(mBest)) + 2 #add 2 parameters to the coefficient count, 1 for variance 1 for delta
  n 		 	<- length(kite.pgr)
  
  mBest$delta 	<- delta.vec
  mBest$LLprof 	<- llVec
  mBest$BIC 	<- -2*logLik(mBest)[1] + k*log(n)
  mBest$AICc  <- -2*logLik(mBest)[1] + 2*k + 2*k*(k+1)/(n-k-1)

  mBest2 	 	<- ll2List[[which.max(ll2Vec)]]
  k 		 	  <- length(coef(mBest2[[1]])) + 2 #add 2 parameters to the coefficient count, 1 for variance 1 for delta
  n 		 	  <- length(kite.pgr)
  
  mBest2$delta 	<- delta.vec
  mBest2$LLprof <- llVec
  mBest2$BIC 	<- -2*logLik(mBest2[[1]])[1] + k*log(n)
  mBest2$AICc  <- -2*logLik(mBest2[[1]])[1] + 2*k + 2*k*(k+1)/(n-k-1)
  
  return(list(mBest, mBest2))

}

mismatch.constrained.optim <- function(kite.pgr, kite.lag, water.cov) {

    delta.vec <- 1/seq(7, 365*2, by=1)
    llVec   <- vector('numeric', length(delta.vec))
    llList  <- list()
    ll2Vec  <- vector('numeric', length(delta.vec))
    ll2List <- list()
  
    for(i in 1:length(delta.vec)) {
      llList[[i]] <- ricker.decay(delta.val=delta.vec[i], kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=water.cov)
      llVec[i]    <- logLik(llList[[i]])[1]

      ll2List[[i]] <- ricker.decaysq(delta.val=delta.vec[i], kite.pgr=kite.pgr, kite.lag=kite.lag, water.cov=water.cov)
      ll2Vec[i]    <- logLik(llList[[i]])[1]

    }

  mBest     <- llList[[which.max(llVec)]]
  k       <- length(coef(mBest)) + 2 #add 2 parameters to the coefficient count, 1 for variance 1 for delta
  n       <- length(kite.pgr)
  
  mBest$delta   <- delta.vec
  mBest$LLprof  <- llVec
  mBest$BIC   <- -2*logLik(mBest)[1] + k*log(n)
  mBest$AICc  <- -2*logLik(mBest)[1] + 2*k + 2*k*(k+1)/(n-k-1)

  mBest2    <- ll2List[[which.max(ll2Vec)]]
  k       <- length(coef(mBest2)) + 2 #add 2 parameters to the coefficient count, 1 for variance 1 for delta
  n       <- length(kite.pgr)
  
  mBest2$delta  <- delta.vec
  mBest2$LLprof <- llVec
  mBest2$BIC  <- -2*logLik(mBest2)[1] + k*log(n)
  mBest2$AICc  <- -2*logLik(mBest)[1] + 2*k + 2*k*(k+1)/(n-k-1)

  return(list(mBest, mBest2))

}

#calculate the ricker model with weighted temperature coveriate, with fixed decay rate
ricker.decay.constraint <- function(par.vec, kite.pgr, kite.lag, water.cov) {
    a <- par.vec[1]
    b <- par.vec[2]
    d <- par.vec[3]
    sig <- exp(par.vec[4])

    delta.val <- (1 - b)/365
    #print(delta.val)
#cat(b, delta.val, '\n')
#stop()
    WeightedCov   <- vector('numeric', length(kite.pgr))
    for(j in 1:length(kite.pgr)) {

      currdays        <- ((j-1)*365+1):(j*365)
      weightVec       <- (1 - delta.val)^(365:1-1)
      WeightedCov[j]  <- sum(weightVec*(water.cov[currdays]))

    }#end for j

    WeightedCov <- (WeightedCov-mean(WeightedCov))/sd(WeightedCov)

    pgr.predict <- a + b*kite.lag + d*WeightedCov
    x <- dnorm(kite.pgr, pgr.predict, sig, log=T)

    return(-sum(x))

}

ricker.decay.constraint.posthoc <- function(par.vec, kite.pgr, kite.lag, water.cov, discrete.cov, water.rmvec) {
    a <- par.vec[1]
    b <- par.vec[2]
    d <- par.vec[3]
    e <- par.vec[4]
    sig <- exp(par.vec[5])
    ind <- 1:length(kite.lag)
    
    if(!is.null(water.rmvec)) {discrete.cov <- discrete.cov[water.rmvec]}
    delta.val <- (1 - b)/365

    WeightedCov   <- vector('numeric', length(kite.pgr))
    for(j in 1:length(kite.pgr)) {

      currdays      <- ((j-1)*365+1):(j*365)
      weightVec       <- (1 - delta.val)^(365:1-1)
      #weightVec       <- exp(-delta.val*(365:1-1))
      WeightedCov[j]  <- sum(weightVec*(water.cov[currdays]))

    }#end for j
    WeightedCov        <- (WeightedCov-mean(WeightedCov))/sd(WeightedCov)
    discrete.cov[ind]  <- (discrete.cov[ind] - mean(discrete.cov[ind]))/sd(discrete.cov[ind])

    pgr.predict <- a + b*kite.lag + d*WeightedCov + e*discrete.cov
    
    x <- dnorm(kite.pgr, pgr.predict, sig, log=T)

    return(-sum(x))

}

ricker.decay.constraintSQ <- function(par.vec, kite.pgr, kite.lag, water.cov) {

    a <- par.vec[1]
    b <- par.vec[2]
    d <- par.vec[3]
    d2 <- par.vec[4]
    sig <- exp(par.vec[5])

    delta.val <- (1 - b)/365

    WeightedCov <- Weighted2Cov  <- vector('numeric', length(kite.pgr))
    for(j in 1:length(kite.pgr)) {

      currdays      <- ((j-1)*365+1):(j*365)
      weightVec       <- (1 - delta.val)^(365:1-1)
      #weightVec       <- exp(-delta.val*(365:1-1))
      WeightedCov[j]  <- sum(weightVec*(water.cov[currdays]))
      Weighted2Cov[j] <- sum(weightVec*(water.cov[currdays]^2))

    }#end for j

    #return(lm(kite.pgr ~ kite.lag + WeightedCov))
    pgr.predict <- a + b*kite.lag + d*WeightedCov + d2*Weighted2Cov
    x <- dnorm(kite.pgr, pgr.predict, sig, log=T)

    return(-sum(x))

}


#calculate the ricker model with weighted temperature coveriate, with fixed decay rate
ricker.decay <- function(delta.val, kite.pgr, kite.lag, water.cov) {

  	WeightedCov 	<- vector('numeric', length(kite.pgr))
	#plot(weightVec, type='l')
  	for(j in 1:length(kite.pgr)) {

    	currdays  		<- ((j-1)*365+1):(j*365)
    	weightVec     	<- (1 - delta.val)^(365:1-1)
    	weightVec 		<- exp(-delta.val*(365:1-1))
    	WeightedCov[j]  <- sum(weightVec*(water.cov[currdays]))

  	}#end for j
  	WeightedCov <- (WeightedCov-mean(WeightedCov))/sd(WeightedCov)

  	return(lm(kite.pgr ~ kite.lag + WeightedCov))

}


#calculate the ricker model with weighted temperature coveriate, with fixed decay rate
ricker.decaysq <- function(delta.val, kite.pgr, kite.lag, water.cov) {

  WeightedCov 	<- vector('numeric', length(kite.pgr))
  WeightedCov2 	<- vector('numeric', length(kite.pgr))

  for(j in 1:length(kite.pgr)) {

    currdays  		  <- ((j-1)*365+1):(j*365)
    weightVec     	<- (1 - delta.val)^(365:1-1)
    WeightedCov[j]  <- sum(weightVec*(water.cov[currdays]))
	  WeightedCov2[j] <- sum(weightVec*(water.cov[currdays]^2))

  }#end for j

  return(list(lm=lm(kite.pgr ~ kite.lag + WeightedCov + WeightedCov2), WeightedCov=WeightedCov, WeightedCov2=WeightedCov2))

}


#analyze the snail kite example
kite.est    <- as.matrix(read.csv(file="../Kite/SNKI abundance_south.csv",header=T))
#kite.est    <- kite.est[1:12,]
years.kite  <- kite.est[,1]

WCA3A.dat   <- read.csv(file="../Kite/WaterData2.csv",header=T,sep='\t')
WCA3A.dat   <- as.numeric(as.vector(WCA3A.dat[,3]))
water.dat   <- WCA3A.dat

##fix the NA's in water level data
days    <- 1:length(WCA3A.dat)
lo 	    <- loess(WCA3A.dat ~ days)
na.vec  <- which(is.na(water.dat))

WCA3A.dat[na.vec] <- lo$fitted[na.vec]
water.dat[na.vec] <- lo$fitted[na.vec]

#get average yearly water level
water.mean		<- vector('numeric', length(kite.est[,1])+1)
water.mean    <- vector('numeric', length(kite.est[,1])+1)
water.rec <- water.max <- water.min <- water.sd <- water.med <- water.dam <- water.mean

for(i in 1:14) {
  ind <- ((i-1)*365+1):(365*i)
  curr.water      <- water.dat[ind]
  water.mean[i]   <- mean(curr.water, na.rm=T)
  water.rec[i]    <- (curr.water[1] -  min(curr.water, na.rm=T))/which.min(curr.water)
  water.max[i]    <- max(curr.water, na.rm=T)
  water.min[i]    <- min(curr.water, na.rm=T)
  water.sd[i]     <- sd(curr.water, na.rm=T)  
  water.med[i]    <- median(curr.water, na.rm=T)  
  water.dam[i]    <- sum(curr.water > water.mean[i])/length(curr.water)
}

kite.lag    <- log(kite.est[1:(length(years.kite)-1),2])
kite.pgr    <- log(kite.est[2:(length(years.kite)),2]/(kite.est[1:(length(years.kite)-1),2]))
rmyear.vec	<- -c(1:2) #remove first and last year from water covariate data to match population growth data

ind         <- 1:10
mean.fit		<- lm.fits(kite.pgr=kite.pgr[ind], kite.lag=kite.lag[ind], water.cov=water.mean, water.rmvec=rmyear.vec)

#rsquared 0. 0.2188201 0.2305357
rec.fit     <- lm.fits(kite.pgr=kite.pgr[ind], kite.lag=kite.lag[ind], water.cov=water.rec, water.rmvec=rmyear.vec)
max.fit     <- lm.fits(kite.pgr=kite.pgr[ind], kite.lag=kite.lag[ind], water.cov=water.max, water.rmvec=rmyear.vec)
min.fit     <- lm.fits(kite.pgr=kite.pgr[ind], kite.lag=kite.lag[ind], water.cov=water.min, water.rmvec=rmyear.vec)
sd.fit      <- lm.fits(kite.pgr=kite.pgr[ind], kite.lag=kite.lag[ind], water.cov=water.sd, water.rmvec=rmyear.vec)
med.fit     <- lm.fits(kite.pgr=kite.pgr[ind], kite.lag=kite.lag[ind], water.cov=water.med, water.rmvec=rmyear.vec)
dam.fit     <- lm.fits(kite.pgr=kite.pgr[ind], kite.lag=kite.lag[ind], water.cov=water.dam, water.rmvec=rmyear.vec)

BIC.vec <- c("median"=unname(med.fit[3]), "mean"=unname(mean.fit[3]), "rec"=unname(rec.fit[3]), "max"=unname(max.fit[3]), "min"=unname(min.fit[3]), "SD"=unname(sd.fit[3]), "DAM"=unname(dam.fit[3]))


#fit the mismatch model
library(nloptr) 
lm.mod <- lm(kite.pgr[ind] ~ kite.lag[ind])#initial parameter estimates
par.vec   <- c(a=coef(lm.mod)[1], b1=coef(lm.mod)[2], d1=0, sigma=log(sd(kite.pgr)))
#fit it
mismatch.constraint   <- mlsl(x0=par.vec, fn=ricker.decay.constraint, lower=c(a=-10, b1=-10, d1=-10, sig=-5), upper=c(a1=10, b1=10, d1=10, sig=1), local.method = "LBFGS", low.discrepancy = TRUE, nl.info = FALSE, control = list(maxeval=1e3), kite.pgr=kite.pgr[ind], kite.lag=kite.lag[ind], water.cov=log(WCA3A.dat[-(1:(182 + 1*365))]))

#now get the standard errors
library(numDeriv)
hmat <- hessian(func=ricker.decay.constraint, x=mismatch.constraint$par, kite.pgr=kite.pgr[ind], kite.lag=kite.lag[ind], water.cov=log(WCA3A.dat[-(1:(182+1*365))]), method.args=list(eps=1e-4, d=1e-6))
semat <- solve(hmat)

names(mismatch.constraint$par) <- names(par.vec)
k <- length(mismatch.constraint$par)
n <- max(ind)
mismatch.constraint$BIC <- 2*mismatch.constraint$value + k*log(n)

# AICc.vec-mismatch.constraint$AICc
#   median      mean       rec       max       min        SD       DAM 
# 7.048186  7.251049  6.208120  6.875634 -1.707940 -1.232502  7.286099

print(mismatch.constraint)

#get Rsquared of mismatch model
water.cov <- log(WCA3A.dat[-(1:(182 + 1*365))])

a <- mismatch.constraint$par[1]
b <- mismatch.constraint$par[2]
d <- mismatch.constraint$par[3]

delta.val <- (1 - b)/365
WeightedCov   <- vector('numeric', length(kite.pgr[1:10]))

for(j in 1:length(kite.pgr[1:10])) {

  currdays        <- ((j-1)*365+1):(j*365)
  weightVec       <- (1 - delta.val)^(365:1-1)
  WeightedCov[j]  <- sum(weightVec*(water.cov[currdays]))

}#end for j

WeightedCov <- (WeightedCov-mean(WeightedCov))/sd(WeightedCov)

pgr.predict <- a + b*kite.lag[1:10] + d*WeightedCov

Rsq <- 1 - var(pgr.predict - kite.pgr[1:10])/var(kite.pgr[1:10])
