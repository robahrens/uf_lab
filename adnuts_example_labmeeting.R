#### MCMC with the adnuts R package
#### requires ADMB version 12.0

## load packages
library(adnuts)
library(coda)
library(rstan)
library(R2admb) # for setting up admb environment, if necessary

as.numeric.factor <- function(x) {as.numeric(as.character(x))}

## set path to folder with admb files (tpl, dat, ctl files)
path <- "/Users/cfriess/Documents/hake_pm"
setwd(path)

## might need to set up admb environment with the R2admb function setup_admb
# setup_admb()

## build executable
system("admb pm")
## run the model
system("./pm")

## run some random walk markov exploratory chains (rwm) and inspect pairs
## want to start these from MLE estimates. can read these in from std file
std.file <- read.table(paste("pm.std",sep=""))
std.file
# get npar - exclude depletion (is not an estimated par) and first row (var names)
lrow.std <- which(std.file[,2] == "depletion") - 1
npar_est <- lrow.std - 1
npar_est
par.values <- as.numeric.factor(std.file[2:lrow.std,3]) 
par.sd <- as.numeric.factor(std.file[2:lrow.std,4])

# run adnuts::sample_admb with rwm option
##########################################

nchains_rwm <- 5 # specify number of chains to run
## inital par values for chains must be in list format
nuts.init_rwm <- rep( list(as.list(par.values)), nchains_rwm ) 
## want to aim for about 1000 samples
iter_rwm <- 10000
thin_rwm <- 10
set.seed(453) # optional to get same random number generator start chains
# if no warmup value is provided, the default is to use half of iter
# can pass extra args to admb as per below. -mcrb reduces par correlation
fit_rwm <- adnuts::sample_admb(model="./pm", init=nuts.init_rwm, 
    path=path, chains = nchains_rwm, seeds = sample(1:1000,nchains_rwm), 
    iter = iter_rwm, thin = thin_rwm, warmup = 2000,
    algorithm = "RWM", extra.args = "-mcrb 5")

### NOTE: can also run chains in parallel using the snowfall package
# we're not doing that here

# inspect results - leading parameters only here
# trace is the default for diag, other options are acf and hist
adnuts::pairs_admb(fit_rwm, diag = "acf", pars = c(1:5)) # autocorrelation
adnuts::pairs_admb(fit_rwm, diag = "trace", pars = c(1:5))
# check for problem signs, e.g. slowly mixing parameters, pars estimated at the bounds
# if necessary, fix model, e.g. by specifying more informative priors,
# or fixing pars that are estimated at bounds, or reparameterizing model



# run adnuts::sample_admb with nuts option (hameltonian mcmc with no-u-turn sampling)
##########################################

# set starting values - want these to be overdispersed
nchains <- 3
nuts.init <- rep( list(list()), nchains ) 
set.seed(98) # optional
for (i in 1:nchains){
    nuts.init[[i]] <- as.list(rnorm(length(par.values), par.values, par.sd*2))
}
# unlist(nuts.init[[1]]) # to inspect 
# if any of these happen to be outside of par bounds, 
# admb automatically puts them just inside bound so no worries

iter = 1000 # default value is 2000
set.seed(346) # optional to get the same seeds for repeated runs
seed_nuts <- sample(1:1000,3)
# adapt_delta of 0.8 is the default - this is the target acceptance rate
# it needs to be increased to get rid of divergent transitions
# the default warmup is half of iter. during this time the algorithm
# adapts the step size
fit_nuts <- sample_admb(model="./pm", init=nuts.init, path=path, 
	chains = nchains, seeds = seed_nuts, iter = iter,
    control = list(adapt_delta = .8))

## look at results!
## can check out specific call(s) to admb with:
fit_nuts$cmd
sp <- extract_sampler_params(fit_nuts, inc_warmup = FALSE)
head(sp) # max treedepth is set to 12. can be increased manually
# treedepth prevents excessively long trajectories, which may be needed
# for poorly specified models
ndiv <- sum(sp$divergent__ == 1) # number of divergent transitions (want 0)
ndiv
pairs_admb(fit_nuts, pars = c(1:5)) # the red dots are divergent transitions
# well, all but one which is the MLE estimate. perhaps poor choice of color)

# rerun with higher target acceptance rate ( = force smaller step sizes)
fit_nuts_d95 <- sample_admb(model="./pm", init=nuts.init, path=path, 
	chains = nchains, seeds = seed_nuts, iter = iter,
    control = list(adapt_delta = .95))

fit <- fit_nuts_d95

pairs_admb(fit, pars = c(1:5))
sp <- extract_sampler_params(fit, inc_warmup = FALSE)
sum(sp$divergent__ == 1) # much improved!

### performance - using rstan package (can also use coda - gives different values)
mon <- rstan::monitor(fit$samples, print = FALSE)
which.max(mon[(1:nrow(mon)-1),'Rhat']) # don't include the log posterior likelihood
minESS <- min(mon[1:npar_est,'n_eff']) # smallest effective sample size
minESS
tot_time <- sum(fit$time.total / 60)
perf <- minESS / tot_time # sampling efficiency performance metric

### plot posteriors and priors
source("read.admb.R")
source("kobe.R")
source("mcpars.R")

LMC=read.psv("pm.psv") # read psv file
# backtransform pars; select only leading pars
LMC=cbind(exp(LMC[,1]),LMC[,2],exp(LMC[,3]),sqrt(1./exp(LMC[,4])),sqrt(1./exp(LMC[,5])))
lbs=c("K","r","q","sig","tau")	
plot.pairs(LMC,label=lbs,prior=c(1,1,0,2,2),mome1=c(8,-1.38,0,3.79,1.71),mome2=c(0.25,0.51,0,0.0102,0.0086))


### can also run mceval on saved samples
system("./pm -mceval -noest -nohess")


# useful info: https://cran.r-project.org/web/packages/adnuts/vignettes/adnuts.html

