# Robert Ahrens
# January 27, 2018
# simple lm 
# 1) simulate the data
# 2) estimate the trend using STAN

# inits
	set.seed(555)
	nobs = 10
	tm=0.2
	tb=4

# simulate data
	x = runif(nobs)*100
	y = tb+tm*x+rnorm(nobs,0,4)

# plot data
	plot(x,y,pch=19,col="steelblue")
###################################################################################

# Estimate using STAN

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# options(mc.cores = 3)

# declare stan settings
	chains = 4
	iter = 1000
	warmup = 0.5*iter # burn-in
	thin = 1 # 1 is recommended
	seed = 777 # default: sample.int(.Machine$integer.max, 1)

# prep stan.data, as named list
	stan.data = list(x = x,
					 y = y,
					 nobs = nobs)

# prep stan.inits, as named list
	inits.func = function()
	{ 
		list(m = rnorm(1,tm, 1), 
			 b = rnorm(1,tb, 3), 
			 sigma = runif(1, min = 0, max = 8))
	}
	stan.inits = replicate(chains,inits.func(),simplify=FALSE)

# run stan
	file.name = c("/Users/robertahrens/Documents/Work/StanExamples/lm.stan") # path to stan file
	fit = stan(file = file.name,
		       data = stan.data,
		       init = stan.inits,
		       chains = chains,
		       warmup = warmup,
		       iter = iter,
		       thin = thin,
		       seed = seed,
		       control = list(adapt_delta = 0.8,max_treedepth=10))

# STAN diagnostics
	print(fit, digits = 3)	
	check_hmc_diagnostics(fit)
	traceplot(fit, pars = c("m","b", "sigma"), inc_warmup = TRUE, nrow = 4)
	pairs(fit, pars = c("m","b", "sigma"))		

# shinystan diagnostic GUI
	# library(shinystan)
	# launch_shinystan(fit)

# extract estimated mu
	mu = extract(fit)$mu
	quants.mu = apply(mu,2,quantile,probs=c(0.25,0.5,0.975))

# plot true with estimate
		require(scales)
		plot(1,1,xlim=c(0,n.years),ylim=c(0,max(quants.mu)),xlab="Years",ylab="Abundance",type="n",lwd=2)
		polygon(c(1:n.years,n.years:1),c(quants.mu[1,],rev(quants.mu[3,])),col=alpha("dodgerblue1",0.4),border="dodgerblue1")
		lines(quants.mu[2,],lwd=3,col="dodgerblue1") # estimate
		lines(N.vec,lwd=3,col="darkorange") # true

