#Mismatched time series model for snail kite
#conducts analysis 
library(nloptr) #does optimization
library(AICcmodavg) #to calculate BIC

#fits all variations of the linear model including density independent, density dependent, and density dependent with covariate
#returns a list that contains the AIC values, log-liklihood values, and BIC values 
#pgr.obs is the observed population growth rate
#lag.obs is the observed lag abundance prediction
#env.cov is the predictor variable
lm.fits <- function(pgr.obs, lag.obs, env.cov) {

  #fits each of the models
  DI.lm     <- lm(pgr.obs ~ 1) #density independence
  DD.lm     <- lm(pgr.obs ~ lag.obs) #denstity dependent
  DDcov.lm  <- lm(pgr.obs ~ lag.obs + env.cov) #density dependent with environmental covariate

  #calculate some summaries of the models
  aicc.table <- c(AIC(DI.lm), AIC(DD.lm), AIC(DDcov.lm)) 
  
  LL.table <- c(logLik(DI.lm), logLik(DD.lm), logLik(DDcov.lm))
  names(LL.table) <- c('Intercept','DD','DD + Mean covariate')

  BIC.table <- c(BIC(DI.lm), BIC(DD.lm), BIC(DDcov.lm))
  names(BIC.table) <- c('Intercept','DD','DD + Mean covariate')
  
  #list of things to return
  return.list         <- list()
  return.list$AICc    <- aicc.table
  return.list$logLik  <- LL.table
  return.list$BIC     <- BIC.table
  
  return(return.list)
}



#function to fit the ricker model with weighted temperature coveriate and with memory parameter that is fixed by the strenght of density dependence as described in Ferguson et al.
#par.vec is a vector of initial parameter values

#a is the intrinsic growth rate, b is the density dependence, d is the impact of the environmental covariate, sig is the (log) standard deviation
#pgr.obs is the observed population growth rate
#lag.obs is the observed lag abundance predictor 
#env.cov is the observed environmental variable
ricker.memory <- function(par.vec, pgr.obs, lag.obs, env.cov) {
    
    a <- par.vec[1]
    b <- par.vec[2]
    d <- par.vec[3]
    sig <- exp(par.vec[4])

    beta.val   <- (1 - b)/365 #this is the memory parameter , delta =1/365
    WeightedCov <- vector('numeric', length(pgr.obs)) #vector of weighted averages of the environmental covariate

    #go through and calculate the 
    for(j in 1:length(pgr.obs)) {
      currdays        <- ((j-1)*365+1):(j*365) #indices for a year of observations
      weightVec       <- (1 - beta.val)^(365:1-1)#weights over the year
      WeightedCov[j]  <- sum(weightVec*(env.cov[currdays])) #weighted average of the environemtn
    }#end for j
    
    WeightedCov <- (WeightedCov - mean(WeightedCov))/sd(WeightedCov) #standardizes the covariate

    predict.pgr <- a + b*lag.obs + d*WeightedCov#calculate predicted pgr

    x <- dnorm(pgr.obs, predict.pgr, sig, log=T) #log-likelihood

    return(-sum(x)) #returns the negative log-likelihood

}

#analyze the snail kite example
kite.est    <- read.csv(file="SNKIabundance.csv",header=T)
N           <- dim(kite.est)[1] #lenght of the time series

WCA3A.dat   <- read.csv(file="WaterData.csv",header=T)
WCA3A.dat   <- log(as.numeric(as.vector(WCA3A.dat[,2])))

days    <- 1:length(WCA3A.dat)
lo 	    <- loess(WCA3A.dat ~ days)
na.vec  <- which(is.na(WCA3A.dat))
WCA3A.dat[na.vec] <- lo$fitted[na.vec]

#form the kite data
lag.obs    <- (kite.est[1:(length(kite.est$year)-1),]$N)
pgr.obs    <- log(kite.est[2:length(kite.est$year),]$N/(kite.est[1:(length(kite.est$year)-1),]$N))
num.years    <- 1:10 #breakpoint occurs post-2007 so only include data up to 2007

#get average yearly water level
water.mean    <- vector('numeric', length(kite.est$year)+1)
for(i in 1:14) {
  ind <- ((i-1)*365+1):(365*i)
  curr.water      <- WCA3A.dat[ind]
  water.mean[i]   <- mean(curr.water, na.rm=T)
}

#fit the density independent model, standard ricker model, and ricker model with covariate. 
mean.fit    <- lm.fits(pgr.obs=pgr.obs[num.years], lag.obs=lag.obs[num.years], env.cov=water.mean[num.years+1])

#guess initial parameters for the mismatch model
lm1.mod   <- lm(pgr.obs[num.years] ~ lag.obs[num.years])
par.vec   <- c(a=coef(lm1.mod)[1], b1=coef(lm1.mod)[2], d1=0, sigma=log(sd(pgr.obs)))

#fit the mismatched timeseries model using nloptr
opts <- list("algorithm"="NLOPT_LN_NELDERMEAD", "xtol_rel"=1.0e-8)
mismatch.fit <- nloptr(x0=par.vec, eval_f=ricker.memory, opts=opts, pgr.obs=pgr.obs[num.years], lag.obs=lag.obs[num.years], env.cov=WCA3A.dat[-(1:(182 + 1*365))])

#AIC and BIC values for the memory model
names(mismatch.fit$solution) <- names(par.vec)
k <- length(mismatch.fit$solution)
n <- max(num.years)

mismatch.fit$BIC <- 2*mismatch.fit$objective + k*log(n)
mismatch.fit$AIC <- 2*mismatch.fit$objective + 2*k 

names(mismatch.fit$BIC) ="DD + Geometric Weighted Covariate"
print(mean.fit$BIC)
print(mismatch.fit$BIC)

