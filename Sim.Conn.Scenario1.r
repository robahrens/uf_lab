

# Nicholas Ducharme-Barth
# January 27, 2018
# simulate multiple noisy indices from: 
# Conn 2010, Hierarchichal analysis of multiple noisy abundance indices
# 1) simulate the data
# 2) estimate the trend using STAN

# inits
	set.seed(999)
	n.years = 20
	n.indices = 10
	q = 0.01

# simulate mean trend
	N.vec = rep(NA,n.years)
	N.vec[1] = 100
	for(i in 2:n.years)
	{
		N.vec[i] = max(c(20,N.vec[i-1] + rnorm(1,0,15^2)))
	}


# simulate noisy indices
# lognormal error structure
	U.mat = matrix(NA,nrow=n.indices,ncol=n.years) # simulated indices
	s.mat = matrix(NA,nrow=n.indices,ncol=n.years) # sample error
	p.vec = rep(NA,n.indices) # process error

	for(i in 1:n.indices)
	{
		U.p = p.vec[i] = sqrt(log(runif(1,0,0.3)^2+1)) # process error
		for(t in 1:n.years)
		{
			U.mean = log(N.vec[t]) + log(q)
			U.s = s.mat[i,t] = sqrt(log(runif(1,0.1,0.3)^2+1)) # sample error

			U.mat[i,t] = exp(rnorm(1,U.mean,U.s+U.p))
		}
	}

# plot noisy data
		plot(1,1,xlim=c(0,n.years),ylim=c(0,max(U.mat)),xlab="Years",ylab="Abundance",type="n",lwd=2)
		for(i in 1:n.indices)
		{
			lines(U.mat[i,],col="gray50")
		}
		lines(q*N.vec,lwd=3,col="darkorange")

###################################################################################
###################################################################################
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
	stan.data = list(U = U.mat,
					 samp_errors = s.mat,
					 n = nrow(U.mat),
					 y = ncol(U.mat))

# prep stan.inits, as named list
	inits.func = function()
	{ 
		list(log_mu = rnorm(ncol(U.mat),log(100), 1), 
			 log_q = rnorm(nrow(U.mat),log(0.01), 0.5), 
			 proc_errors = runif(nrow(U.mat), min = 0, max = 5))
	}
	stan.inits = replicate(chains,inits.func(),simplify=FALSE)

# run stan
	file.name = c("/Users/robertahrens/Documents/Work/StanExamples/Conn.scenario.1.stan") # path to stan file
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
	traceplot(fit, pars = c("mu[1]","mu[2]", "q[1]","q[2]", "proc_errors[1]","proc_errors[2]"), inc_warmup = TRUE, nrow = 4)
	pairs(fit, pars = c("mu[1]","mu[2]", "q[1]","q[2]", "proc_errors[1]","proc_errors[2]"))		

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

