# Code and data for modeling mismatched time series

The files in [Writing/ExampleCode](https://github.com/troutinthemilk/WithinBetween/tree/master/Writing/ExampleCode) are the most useful examples to look at and are documented below. Other files in this repository were used in the manuscript [Detecting population–environmental interactions with mismatched time series data](http://onlinelibrary.wiley.com/doi/10.1002/ecy.1966/abstract). 

*GompertzSim.R* - Simulates from a continuous-time Gompertz model driven by a white noise covariate. It then fits several discrete time models with and without covariates. The covariates fit are the standard mean and the geometrically weighted mean. Main functions are Gomp.fit, which profiles over a range of decay values to fit the discrete Gompertz with geometrically weighted covariates, Gomp.decay, which is called by Gomp.fit and returns the log-likelihood for the Gompertz with a geometrically weighted covariate. The last function is gomp.onestep, which simulates from a continuous time Gompertz model using Gauss’s method.

*SnailKiteMismatch.R* - Applies the geometrically weighted covariate model to the snail kite dataset. The function lm.fits fits all variations of the standard linear model with and without covariate to the data. Outputs a list object that contains AIC values, BIC values, and log-likelihood values. The function ricker.memory fits the geometrically weighted covariate model and outputs the negative log-likelihood.

*SNKIabundance.csv* - Data columns are year, the sampling year, and N, the estimated abundance of the southern snail kite population. 


*WaterData.csv* - Data columns are the Daily Date, in Day-Month-Year format, and the Daily Value, which is the daily water level used as a covariate to predict snail kite population growth. 
