
library('MuMIn')
library(abind)
library(foreach)
library(doMC)
options(na.action = "na.fail")


batchFlag   = !is.na(Sys.getenv("R_BATCH", NA))
if(batchFlag) {registerDoMC(4)} else {registerDoMC(1)}

insectData    = read.csv(file='../Data/Deutsch_data.csv')
tempData      = read.csv(file='../Data/TempData.csv')
mintol        = 1e-4
lifetime.vec  = exp(seq(log(7), log(2*365), length.out=64))  
delta.vec     = 1/lifetime.vec

logistic.PT.rmax = function(PT, rmax=1, delta=1, span=1, N0=1) {

  Nvec  	  = vector('numeric', length(PT))
  Nvec[1] 	= N0
  dt  		  = span #/length(PT)

  for(i in 2:length(PT)) {
    Nvec[i] = max(Nvec[i-1] + dt*(rmax*min(max(PT[i-1],0), 1)*Nvec[i-1] - delta*Nvec[i-1]^2), 0)
  }

  return(Nvec)

}

PT.func = function(Temp, Topt, sigmap, CTmax) {
  PT = vector('numeric', length(Temp))
  for(i in 1:length(Temp)) {
    if(Temp[i] < Topt) {
      PT[i] = exp(-((Temp[i] - Topt)/(2*sigmap))^2)
    } else {
      PT[i] = 1-((Temp[i] - Topt)/(Topt - CTmax))^2
    }
  }
  if(any(PT < 0)) {PT[PT<0] = 0}
  return(PT)
}

                       #WeightedT, localT, WeightedT2, avgT
Ricker.fit = function(x, WeightedT, localTemp, WeightedT2, avgTemp) {

  if(any(is.na(x))) { return(rep(-99, 25)) }

  xcov  	= log(x[-length(x)])
  pgr     = log(x[-1]/x[-length(x)])
  pgr     = scale(pgr)[,1]
  xcov    = scale(xcov)[,1]

  mnull = lm(pgr ~ 1)
  m0  = lm(pgr ~ xcov)

  delta.vec = 1/seq(7, 365*2, by=1)
  llVec = vector('numeric', length(delta.vec))
  llList = list()

  for(i in 1:length(delta.vec)) {
    llList[[i]] = ricker.decay(delta.vec[i], pgr=pgr, xcov=xcov)
    llVec[i]    = logLik(llList[[i]])[1]
  }
  m1 = llList[[which.max(llVec)]]
  m1.BIC = -2*logLik(m1)[1] + (length(coef(m1))+2)*log(length(pgr))

  m2  = lm(pgr ~ xcov*WeightedT)
  llVec = vector('numeric', length(delta.vec))
  llList = list()
  for(i in 1:length(delta.vec)) {
    llList[[i]] = ricker.decay2(delta.vec[i], pgr=pgr, xcov=xcov)
    llVec[i] = logLik(llList[[i]])[1]
  }
  m3 = llList[[which.max(llVec)]]
  m3.BIC = -2*logLik(m3)[1] + (length(coef(m1))+3)*log(length(pgr))

  m4  = lm(pgr ~ xcov + WeightedT2) #WeightedT2
  m5  = lm(pgr ~ xcov * WeightedT2)

  llVec = vector('numeric', length(delta.vec))
  llList = list()
  for(i in 1:length(delta.vec)) {
    llList[[i]] = ricker.window(delta.vec[i], pgr=pgr, xcov=xcov)
    llVec[i] = logLik(llList[[i]])[1]
  }
  m6 = llList[[which.max(llVec)]]
  m6.BIC = -2*logLik(m6)[1] + (length(coef(m6))+2)*log(length(pgr))

  m7  = lm(pgr ~ xcov + avgTemp)

  m8  = lm(pgr ~ xcov + WeightedT + WeightedT2 + localTemp)
  
  llVec = vector('numeric', length(delta.vec))
  llList = list()
  for(i in 1:length(delta.vec)) {
    llList[[i]] = rickerniche.decay(delta.vec[i], pgr=pgr, xcov=xcov)
    llVec[i] = logLik(llList[[i]])[1]
  }
  m9 = llList[[which.max(llVec)]]
  m9.BIC = -2*logLik(m9)[1] + (length(coef(m9))+2)*log(length(pgr))

  Rsq.vec <- c(rsq.null=summary(mnull)$r.squared, rsq.m0=summary(m0)$r.squared, rsq.m1=summary(m1)$r.squared, rsq.m2=summary(m2)$r.squared, rsq.m3=summary(m3)$r.squared, rsq.m4=summary(m4)$r.squared, rsq.m5=summary(m5)$r.squared, rsq.m6=summary(m6)$r.squared, rsq.m7=summary(m7)$r.squared, rsq.m8=summary(m8)$r.squared, rsq.m9=summary(m9)$r.squared)
  BIC.vec <- c(BIC.null=BIC(mnull), BIC.m0=BIC(m0), BIC.m1=m1.BIC, BIC.m2=BIC(m2), BIC.m3=m3.BIC, BIC.m4=BIC(m4), BIC.m5=BIC(m5), BIC.m6=m6.BIC, BIC.m7=BIC(m7), BIC.m8=BIC(m8), BIC.m9=BIC(m9))

  return(c(Rsq.vec, BIC.vec, coef(m1)))

}

ricker.decay <- function(parvec, pgr, xcov) {

  delta=parvec[1]

  weightVec = vector('numeric', length(pgr))
  WeightedT = vector('numeric', nlevels(yearfact)-1)

  for(j in 1:(nlevels(yearfact)-1)) {

    currdays	= get.currdays(j, months, yearfact, endMonth, startMonth, indices)

    if(endDayMax > length(currdays)) { endDay = startDay + length(currdays)-1 } else{endDay = endDayMax}

    weightVec     = (1 - delta)^(endDay:startDay-startDay)
    WeightedT[j]   = mean(weightVec*(currtemp[currdays]))

  }#end for j

  return(lm(pgr ~ xcov + WeightedT))

}

rickerniche.decay <- function(parvec, pgr, xcov) {

  delta=parvec[1]

  weightVec = vector('numeric', length(pgr))
  WeightedT = vector('numeric', nlevels(yearfact)-1)
  WeightedT2 = WeightedT

  for(j in 1:(nlevels(yearfact)-1)) {

    currdays  = get.currdays(j, months, yearfact, endMonth, startMonth, indices)

    if(endDayMax > length(currdays)) { endDay = startDay + length(currdays)-1 } else{endDay = endDayMax}

    weightVec      = (1 - delta)^(endDay:startDay-startDay)
    WeightedT[j]   = mean(weightVec*(currtemp[currdays]))
    WeightedT2[j]  = mean(weightVec*(currtemp[currdays]^2))

  }#end for j

  return(lm(pgr ~ xcov + WeightedT + WeightedT2))

}

ricker.decay2 <- function(parvec, pgr, xcov) {

  delta=parvec[1]

  weightVec = vector('numeric', length(pgr))
  WeightedT = vector('numeric', nlevels(yearfact)-1)

  for(j in 1:(nlevels(yearfact)-1)) {

    currdays	= get.currdays(j, months, yearfact, endMonth, startMonth, indices)

    if(endDayMax > length(currdays)) { endDay = startDay + length(currdays)-1 } else{endDay = endDayMax}

    weightVec     = (1 - delta)^(endDay:startDay-startDay)
    WeightedT[j]   = mean(weightVec*(currtemp[currdays]))

  }#end for j

  return(lm(pgr ~ xcov  + WeightedT + I(WeightedT^2)))

}


ricker.window <- function(parvec, pgr, xcov) {

  delta=1/parvec[1]

  weightVec = vector('numeric', length(pgr))
  localT = vector('numeric', nlevels(yearfact)-1)

  for(j in 1:(nlevels(yearfact)-1)) {

    currdays	= get.currdays(j, months, yearfact, endMonth, startMonth, indices)

    if(endDayMax > length(currdays)) { endDay = startDay + length(currdays)-1 } else{endDay = endDayMax}

    dayInd   		   = floor(max(length(currdays)-lifetime.vec[i], 1)):min(length(currdays), 365)
    localT[j]      = mean((currtemp[currdays[dayInd]]))

  }#end for j

  return(lm(pgr ~ xcov + localT))

}


loess.fitpredict <- function(PT.val, currdays, spandays=7, scaleFactor=1) {

  PT.loess 		<- loess(PT.val[currdays]~currdays, span=spandays/length(currdays))
  newdat 		<- seq(currdays[1],currdays[length(currdays)], length.out=length(currdays)*scaleFactor)
  newdat 		<- data.frame(currdays=newdat)

  PT.predicted	<- predict(PT.loess, newdat)

  return(PT.predicted)

}


get.currdays <- function(j, months, years, endMonth, startMonth, indices) {

      yearfact = as.factor(years)
      curryear = levels(yearfact)[j]

      if(endMonth > startMonth) {
        initInd = which(yearfact==curryear & months == startMonth)[1]
        endInd  = which(yearfact==curryear & months == endMonth)
        endInd    = endInd[length(endInd)]
        currdays    = initInd:endInd
      } else {
        nextyear  = levels(yearfact)[j+1]
        initInd = which(years == curryear & months == startMonth)[1]
        endInd  = which(years == nextyear & months == endMonth)
        endInd  = endInd[length(endInd)]
        currdays  = initInd:endInd
      }
      if(length(indices)<366) {
        currdays    = currdays[indices]
        if(any(is.na(currdays))) {currdays = currdays[-which(is.na(currdays))]}
      } else {currdays = currdays[1:length(currdays)]}
      return(currdays)
}


search.func <- function(r.vec, PT.val, months, years, endMonth, startMonth, startDay, endDayMax, indices, currdelta, opt=T) {

  yearfact  = as.factor(years)
  xpt       = vector('numeric', nlevels(yearfact)-1)#matrix(NA, length(r.vec), nlevels(yearfact)-1)

  for(j in 1:(nlevels(yearfact)-1)) {
    currdays   = get.currdays(j, months, years, endMonth, startMonth, indices)
    if(endDayMax > length(currdays)) { endDay = startDay + length(currdays)-1 } else{ endDay = endDayMax }
      PT.predicted  = PT.val[currdays]
      #PT.predicted    = loess.fitpredict(PT.val, currdays, spandays=7, scaleFactor=50)

      if(j==1) { N0=1 } else{ N0=xpt[j-1] }
        temp      = logistic.PT.rmax(PT.predicted, r.vec, delta=currdelta, span=1, N0=N0)
        xpt[j]    = temp[length(temp)]
      } #end for j
    if(opt==TRUE) { return(abs(xpt[length(xpt)]-1)) }
    else { return(xpt) }
}

PT.limits <- function(PT.val, currdays) {

	PT.predicted  	= loess.fitpredict(PT.val, 1:length(PT.val), spandays=30)
	PTmax 			= which.max(PT.predicted)
	tempHi	 	  	= which(PT.predicted <= 0.1)
	startDay		= which(tempHi < PTmax)
	endDay			= which(PT.predicted <= 0.1 & as.numeric(names(PT.predicted)) > PTmax)
	startDay 		= max(startDay[length(startDay)]+1, 1)
	endDay 			= min(endDay, 365, na.rm=T)

	return(list(PT.predicted, startDay, endDay))

}


acomb <- function(...) { abind(..., along=3) }

Rvec      = vector('numeric', length(delta.vec))
Rpgrvec   = Rvec
avgTvec   = Rvec
rmaxbest  = Rvec
meanPtvec = Rvec
varPtvec  = Rvec
random    = FALSE

Rsq.matrix  = foreach(i = 1:length(delta.vec), .combine='acomb', .multicombine=TRUE) %do% {

   cat("i",i,'\n')
   ricker.fits <- matrix(NA, dim(insectData)[1], 25)

   for(h in 1:(dim(insectData)[1])) {

      cat('h', h, '\n')

      if(h == 33) {
  	  cat(' end h', h, '\n')
      ricker.fits[h,] = rep(-99, 25)
      next()
    }

    currtemp  = tempData[,h+1]
    names(currtemp) = tempData[,1]

    sigmap  = (insectData[h,]$Topt - insectData[h,]$Ctmin)/4
    PT.val  = PT.func(currtemp - 273.15, insectData[h,]$Topt, sigmap, insectData[h,]$CTmax)
    Tdates  = as.Date(names(currtemp), "%Y_%m_%d")
    years   = as.numeric(format(Tdates, "%Y"))
    months  = as.numeric(format(Tdates, "%m"))
    yearfact= as.factor(years)

    startMonth  = insectData[h,]$StartMonth
    endMonth    = (insectData[h,]$StartMonth + (11))%%12
    indices     = insectData[h,]$StartDay:insectData[h,]$EndDay
    if(random) { indices = sample(indices) }
    if(endMonth==0) { endMonth=12 }

    ##now resolve to be a finer search until tolerance is reached
    startDay		= insectData[h,]$StartDay
  	endDayMax 	= insectData[h,]$EndDay

  	f.min = optimize(f=search.func, interval=c(0,1), PT.val=PT.val, months=months, years=years, endMonth=endMonth, startDay=startDay, endDayMax=endDayMax, startMonth=startMonth, indices=indices, currdelta=delta.vec[i])
    x     = search.func(r.vec=f.min$minimum, PT.val=PT.val, months=months, years=years, endMonth=endMonth, startMonth=startMonth, startDay=startDay, endDay=endDayMax, indices=indices, currdelta=delta.vec[i], opt=FALSE)

    #x = search.func(r.vec=delta.vec[i], PT.val=PT.val, months=months, years=years, endMonth=endMonth, startMonth=startMonth, startDay=startDay, endDay=endDayMax, indices=indices, currdelta=delta.vec[i], opt=FALSE)

    if(f.min$objective > 0.9) {
		    cat(' end h', h, '\n')
    	cat('failed to converge', i, h, f.min$objective, '\n')
      print(f.min)
    	ricker.fits[h,] <- rep(-99, 25)
    	next()
  	}

  	x	= search.func(r.vec=f.min$minimum, PT.val=PT.val, months=months, years=years, endMonth=endMonth, startMonth=startMonth, startDay=startDay, endDay=endDayMax, indices=indices, currdelta=delta.vec[i], opt=FALSE)

    #get yearly average T
    avgT 	= vector('numeric', nlevels(yearfact)-1)
    Phase	= vector('numeric', nlevels(yearfact)-1)
    Freq	= avgT
    localT  = avgT
    WeightedT = avgT
    WeightedT2 = avgT

    for(j in 1:(nlevels(yearfact)-1)) {

      currdays	= get.currdays(j, months, yearfact, endMonth, startMonth, indices)

  	  if(endDayMax > length(currdays)) { endDay = startDay + length(currdays)-1 } else{endDay = endDayMax}

      #avgT[j]     	= mean(currtemp[currdays[startDay:endDay]])

      weightVec     = (1 - delta.vec[i])^(endDay:startDay-startDay)
      weightVec2		= exp(-delta.vec[i]*seq(1, 0, length.out=endDay-startDay+1))
      avgT[j]       = mean(currtemp[currdays])

      WeightedT[j]   = mean(weightVec*(currtemp[currdays]))
      WeightedT2[j]  = mean(weightVec2*(currtemp[currdays]))
      dayInd   		   = floor(max(length(currdays)-lifetime.vec[i], 1)):min(length(currdays), 365)
      localT[j]      = mean((currtemp[currdays[dayInd]]))

    }#end for j
	  ricker.fits[h,] <- Ricker.fit(c(1,x), WeightedT=WeightedT, localTemp=localT, WeightedT2=WeightedT2, avgTemp=avgT)
    print(ricker.fits[h,])
  cat(' end h', h, '\n')

  } #species loop

  cat('end i', i, '\n')
  print(dim(ricker.fits))
  print(ricker.fits)

  if(any(is.na(ricker.fits))) { ricker.fits[is.na(ricker.fits)] <- -99 }

  return(ricker.fits)

}#delta loop

Rsq.matrix[which(Rsq.matrix==-99)] <- NA
dimnames(Rsq.matrix)[[2]] = c('rsq.null', 'rsq.m0', 'rsq.m1', 'rsq.m2', 'rsq.m3', 'rsq.m4', 'rsq.m5', 'rsq.m6', 'rsq.m7', 'rsq.m8', 'rsq.m9', 'BIC.null', 'BIC.m0', 'BIC.m1', 'BIC.m2', 'BIC.m3', 'BIC.m4', 'BIC.m5', 'BIC.m6', 'BIC.m7', 'BIC.m8', 'BIC.m9', 'b0', 'b1', 'b2')

if(random) {save.image(file="InsectDeathRatesMRRandom_joint.Rdata")} else {
  save.image(file="InsectDeathRatesMR_joint.Rdata")
}

if(batchFlag) { q('no') }
