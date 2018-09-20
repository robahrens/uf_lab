
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#			Advection-Diffusion-Reaction Modeling
#		
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Prepared by Rob Ahrens, Fabio Caltabellotta, and Zach Siders

###-----------------------------------------------------
# Two functions:
# 1) fit_ADR simulates data using the params argument for the diffusion (Dx, Dy) and advection parameters (Ux, Uy) then solves for those parameters using optim.
	# It uses two multinomials to simulate data, 1) draws a random day that the tags were recovered on; and 2) draws a cell using probabilities generated from the ADR model using the given parameters. All tags were released in the center of the grid.
	# From the release point, the change in fish location per unit time is determined by an internal "sim" function to generate the derivatives.
	# These derivatives are passed to ode.2D function (deSolve package) to solve those derivatives at given time steps (seq(0, days_at_large)).
	# The derivative solution is converted to a probability and then used in the second multinomial (probability of a tag being recovered at particular cell on a particular day).

	# A second internal function "poissonlike" uses another internal "sim" function to find the derivatives for each new parameter estimate (Dx, Dy, Ux, and Uy) made by optim. These derivatives are solved again by ode.2D and then the probability of getting a tag in a cell on a given day is calculated. These probabilities are then used in a Poisson likelihood (Hilborn Sparse Multinomial - Hilborn 1990 "Determination of fish movement patterns from tag recoveries using maximum likelihood estimators". CJFAS). The data for this likelihood is the observed returns (simulated previously internally).

# 2) vizer visualizes the output of fit_ADR. To do this, it calculates the derivatives for the parameter estimates from fit_ADR, solves them with ode.2D, then generates the probability surface for the given day in question (argument "day"). 

# **** there is no reaction component in this simulation which encompasses the various mortality components that might occur post-tagging (tag mortality, natural mortality, fishing mortality, etc.)

# depending on your machine you may run into an issue where the lrw argument in the ode.2D function kicks an error. If it does, you just need to make the value of that argument larger. You will likely need to change this in all occurences of ode.2D (there are 3: 1) simulation in fit_ADR, 2) likelihood in fit_ADR, and 3) vizer).  Example: DLSODES- RWORK length is insufficient to proceed. set argument lrw larger than LENRW (=I1), is now: LRW (=I2)

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Packages
library(raster)
library(ReacTran)
library(viridis)
library(mvtnorm)
library(MASS)

####ReacTran Diffusion Model####

# Arguments
ntags <- 100 #how many tags to put out
days_at_large <- 365 #how many days were tagged fish at large
nxcell <- 100 #number of x cells
nycell <- 100 #number of y cells
params <- c(Dx=10, Dy=7, Ux=0.4, Uy=0.1)

#fits to the data using optim
fit_ADR <- function(ntags=100, days_at_large=365, nxcell=100, nycell=100, params){
	
	inv_logit <- function(x){
		exp(x)/(exp(x)+1)
	}
		
	######## SIMULATION
		#### SPATIAL JUNK
			ncell <- nxcell*nycell #number of cells in matrix
			extent <- c(1,nxcell,1,nycell) #extent of matrix
			#raster grid to run over
			r <- raster(extent(extent), nrows=nycell, ncol=nxcell)
			dx <- res(r)[1]
			dy <- res(r)[2]

		#### TEMPORAL
			#which day the returns come in
			day_returns <- rmultinom(1, ntags, prob=inv_logit(rnorm(length(2:days_at_large))))
			#this is actual day they come in
			check_days <- which(day_returns>0)

		#### MARKING
			#dumps all tags in the center
			mark_cells <- matrix(0, nxcell, nycell)
			mark_cells[(ncell/2)-(nxcell/2)] <- ntags
			#marks by cell
			omarks <- r
			omarks[] <- 0

			arr_ind <- which(mark_cells > 0, arr.ind=TRUE)
			cells_w_marks <- cellFromXY(omarks, arr_ind)
			
			#fills the cells that got a fish tagged in them
			for(i in 1:length(cells_w_marks)) omarks[cells_w_marks[i]] <- mark_cells[arr_ind[i,1],arr_ind[i,2]]

		#### RETURNS
			##Define initial conditions
			#t is time steps
			#Init is the matrix Nx * Ny that with #'s of tagged fishies in them
			#params are our diffusion and advection values
			sim<-function(t, Init, params)
			{	
				FISH <- matrix(Init, nxcell, nycell)
				#Transport
				#tran.2D returns the derivatives (change in fish location per unit time)
				Dx <- params[1]
				Dy <- params[2]
				Ux <- params[3]
				Uy <- params[4]

				abl <- 0
				r_grow <- 0

				Tran <- tran.2D(FISH, D.x=Dx, D.y=Dy, dx=dx, dy=dy, v.x=Ux, v.y=Uy, a.bl.x.up=abl, a.bl.x.down=abl, a.bl.y.up=abl, a.bl.y.down=abl)

				#Transport + reaction
				dFISH <- (Tran$dC + r_grow*FISH)

				return(list(dFISH))
			}

			#solves for the diffusion and advection parameters
			out <- ode.2D(y=as.vector(omarks), times=seq(0,days_at_large, by=1), func=sim, parms=params, method="lsodes", lrw=1000000, dimens=c(nxcell,nycell))

			#drops weird negatives
			out[out < 0] <- 0
			#probability of return in a cell at a given day
			simp <- out[2:nrow(out),2:ncol(out)]/rowSums(out[2:nrow(out),2:ncol(out)])

			#storage matrix ndays by ncells
			oreturns <- matrix(0, nrow=length(1:days_at_large), ncol=nxcell*nycell)

			#generates the actual day by cell returns
			for(i in 1:length(check_days)){
				oreturns[check_days[i],] <- rmultinom(1, day_returns[check_days[i]], prob=simp[check_days[i],])
			}


			### the rest just goes to the vizer function
				#which cells in the overall map had returns in them
				cells_w_returns <- unlist(apply(oreturns[check_days,], 1, function(x) {which(x > 0)}))

				return_matrix <- matrix(0, nxcell, nycell)

				for(i in 1:length(cells_w_returns)){
					return_matrix[cells_w_returns[i]] <- return_matrix[cells_w_returns[i]] + 1
				}

				return_coords <- which(return_matrix > 0, arr.ind=TRUE)
			
	######## LIKELIHOOD
	
		tD <- 8
		tUx <- 0.3
		tUy <- 0.05
		#initial diffusion and advection parameter guesses
		theta<-c(log(tD),log(tD), tUx, tUy)
		poissonlike<-function(theta)
		{
			sim2 <-function(t, Init, params)
			{	
				FISH <- matrix(Init, nxcell, nycell)
				#Transport
				#tran.2D returns the derivatives (change in fish location per unit time)
				Dx <- exp(params[1])
				Dy <- exp(params[2])
				Ux <- params[3]
				Uy <- params[4]

				abl <- 0
				r_grow <- 0

				Tran <- tran.2D(FISH, D.x=Dx, D.y=Dy, dx=dx, dy=dy, v.x=Ux, v.y=Uy, a.bl.x.up=abl, a.bl.x.down=abl, a.bl.y.up=abl, a.bl.y.down=abl)

				#Transport + reaction
				dFISH <- (Tran$dC + r_grow*FISH)

				return(list(dFISH))
			}
			#solves for the diffusion and advection parameters
			out2 <- ode.2D(y=as.vector(omarks), times=seq(0,days_at_large, by=1), func=sim2, parms=theta, method="lsodes", lrw=10000000, dimens=c(nxcell,nycell))

			out2[out2 < 0] <- 0
			#probability of return in a cell at a given day
			simp <- out2[,2:ncol(out2)]/rowSums(out2[,2:ncol(out2)])


			NLL <- vector(length=length(check_days))
			for(j in 1:length(check_days))
			{
				NLL[j] <- -1*sum(dpois(oreturns[check_days[j],], simp[check_days[j],], log=TRUE)) #Hilborn sparse multinomial
			}
			#Procedure section/objective
			return(sum(ifelse(NLL==Inf,999,NLL)))
		}
		fit <- optim(theta, poissonlike, method="Nelder-Mead", hessian=FALSE)
		results <- c(Dx=exp(fit$par[1]), Dy=exp(fit$par[2]), Ux=fit$par[3], Uy=fit$par[4])


	package <- list(params=results, returns=return_coords, marks=arr_ind, mark_ras=omarks)

	return(package)
}

fit1 <- fit_ADR(ntags=100, days_at_large=365, nxcell=100, nycell=100, params=params)

#visualizes the probability surface over X day
vizer <- function(fit, day=365, plot=TRUE){
	r <- fit$mark_ras
	res <- res(r)
	dx <- res(r)[1]
	dy <- res(r)[2]
	params <- fit$params
	marks <- fit$marks
	returns <- fit$returns


	sim <-function(t, Init, params)
	{	
		abl = 0
		FISH <- matrix(Init, nxcell, nycell)
		#Transport
		Tran <- tran.2D(FISH, D.x=params[1], D.y=params[2], dx=dx, dy=dy, v.x=params[3], v.y=params[4], a.bl.x.up=abl, a.bl.x.down=abl, a.bl.y.up=abl, a.bl.y.down=abl)

		#Transport + reaction
		dFISH <- (Tran$dC + FISH)

		return(list(dFISH))
	}
	out <- ode.2D(func=sim, y=as.vector(fit$mark_ras), times=c(0,day), parms=params, method="lsodes", lrw=1000000, dimens=c(nxcell,nycell))

	out[out < 0] <- 0
	#probability of return in a cell at a given day
	simp <- out[,2:ncol(out)]/rowSums(out[,2:ncol(out)])

	out.mat <- r
	out.mat[] <- simp[2,]

	if(plot==TRUE){
		# dev.new(width=7.2, height=3.9)
		par(mar=c(3,4,1,1))
		plot(out.mat, col=viridis(200), las=1)
		points(marks[,1], marks[,2], pch="X", cex=1.3) #marks
		arrows(x0 =  marks[,1], y0 = marks[,2], x1 = returns[,1], y1 = returns[,2], angle=25, length=0.05)
	}
	return(out.mat)
}
viz <- vizer(fit=fit1, day=30)
