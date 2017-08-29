#simulate from a continuous-time Gompertz model driven by a white noise covariate. 
#fit discrete time model and compare different covariate values

alpha.vec <- 365/exp(seq(log(7), log(2*365), length.out=64)) #vector of timescales to examine. reduce length.out to 


#fit the discrete Gompertz model to discretely sampled data
Gomp.fit <- function(x, env.state) {

  if(any(is.na(x))) { return(rep(-99, 25)) }

  xcov    <- log(x[-length(x)])
  pgr     <- log(x[-1]/x[-length(x)])
  
  delta.vec <- 1/seq(7, 2*365, by=1)
  llVec 	<- vector('numeric', length(delta.vec))
  llList 	<- list()

  for(i in 1:length(delta.vec)) {
    llList[[i]] <- Gomp.decay(delta.vec[i], env.state=env.state, pgr=pgr, xcov=xcov)
    llVec[i]    <- logLik(llList[[i]])[1]
  }

  m1 		<- llList[[which.max(llVec)]]
  m1$BIC 	<- -2*logLik(m1)[1] + (length(coef(m1))+2)*log(length(pgr))

  return(m1)
}

#calculate the gompertz model with weighted temperature coveriate, with fixed decay rate
Gomp.decay <- function(parvec, env.state, pgr, xcov) {

  delta 	<- parvec[1]
  WeightedT <- vector('numeric', length(xcov))

  for(j in 1:length(xcov)) {

    currdays  	  <- ((j-1)*364+1):(j*364)
    weightVec     <- (1 - delta)^(364:1-1)
    WeightedT[j]  <- mean(weightVec*(env.state[currdays]))

  } #end for j

  return(lm(pgr ~ xcov + WeightedT))

}

#simulate from the continous time model using the delta method. projects one step forward
gomp.onestep <- function(alpha, x0=1, steps=365, env.state, delta=1/365) {
	x.vec 		<- vector('numeric', steps)
	x.vec[1] 	<- x0

	for(i in 2:length(x.vec)) {
		x.vec[i] <- x.vec[i-1] + delta*(alpha - alpha/1*x.vec[i-1] + env.state[i-1])
	}

	return(x.vec)
}



##run everything.

#vectors to store rsq values
rsq.gomp 	<- vector('numeric', length(alpha.vec))
rsq.decay 	<- rsq.gomp
rsq.mean	<- rsq.gomp

for(i in 1:length(alpha.vec)) {
	set.seed(1)

	ts.vec 		<- vector('numeric', 1e3)
	ts.vec[1] 	<- 1
	env.state	<- NULL
	env.mean	<- ts.vec

	for(j in 2:length(ts.vec)) {	
		curr.state 	<- rnorm(364, 0, 1)
		env.state	<- c(env.state, curr.state)
		env.mean[j] <- mean(curr.state)
		x 			<- gomp.onestep(alpha=alpha.vec[i], x0=ts.vec[j-1], env.state=curr.state)
		ts.vec[j] 	<- x[length(x)]
	}

	pop.state <- exp(ts.vec)
	
	#fit models 
	pgr 	<- log(pop.state[-1]/pop.state[-length(pop.state)])
	xcov 	<- log(pop.state[-length(pop.state)])
	x 		<- Gomp.fit(x=pop.state, env.state=env.state)

	#store rsqaured values
	rsq.decay[i] <- summary(x)$r.squared
	rsq.gomp[i]  <- summary(lm(pgr~xcov))$r.squared
	rsq.mean[i]  <- summary(lm(pgr~xcov + env.mean[-1]))$r.squared

	cat(i, round(rsq.gomp[i], 2), round(rsq.mean[i], 2), round(rsq.decay[i], 2), '\n')

}


#plot the output
plot(365/alpha.vec, rsq.gomp, ylim=c(0,1), type='l', lwd=3, log="x",xlab="Reproductive period [days]", ylab=expression(paste("Goodness of fit ", (R^2))), col='firebrick3', cex.lab=1.3)
lines(365/alpha.vec, rsq.mean, col="royalblue2", lwd=3)
lines(365/alpha.vec, rsq.decay, col="chartreuse4", lwd=3)
legend('bottomleft', legend=c("Density dependence only", "Mean Environmental covariate", "Weighted environmental covariate"), lwd=3, col=c("firebrick3","royalblue2","chartreuse4"))



#plot the output
#BW vsersion
plot(365/alpha.vec, rsq.gomp, ylim=c(0,1), type='l', lwd=3, log="x",xlab="Reproductive period [days]", ylab=expression(paste("Goodness of fit ", (R^2))), col=rgb(200/256,200/256,200/256), cex.lab=1.3)
lines(365/alpha.vec, rsq.mean, col=rgb(115/256,115/256,115/256), lwd=3)
lines(365/alpha.vec, rsq.decay, col=rgb(0,0,0), lwd=3)
legend('bottomleft', legend=c("Density dependence only", "Mean Environmental covariate", "Weighted environmental covariate"), lwd=3, col=c(rgb(200/256,200/256,200/256), rgb(115/256,115/256,115/256), 'black'))

#library(plotly)
#p <- plot_ly(x=365/alpha.vec, y=rsq.gomp) %>%
#add_trace(x=365/alpha.vec, y=rsq.mean) %>%
#add_trace(x=365/alpha.vec, y=rsq.decay)

# plotly_POST publishes the figure to your plotly account on the web
#plotly_POST(p, filename = "GompertzFits", sharing='public')