
# Nicholas Ducharme-Barth
# February 8, 2018
# EDM (Empirical Dynamic Modeling) Tutorial
# Code modified from rEDM vignette: "rEDM: an R package for Empirical Dynamic Modeling and Convergent Cross Mapping"
# by Hao Ye, Adam Clark, Ethan Deyle, and George Sugihara
# For a more in-depth tutorial please see the many vignettes available at:
# https://cran.r-project.org/web/packages/rEDM/index.html

# install packages (if necessary)
	# install.packages("rEDM")
	# install.packages("scales")
	# install.packages("rgl")

# load packages
	library(rEDM)
	# library(rgl)

# load example data
	data(sardine_anchovy_sst)
	anchovy = sardine_anchovy_sst$anchovy
	sardine = sardine_anchovy_sst$sardine
	sio_sst = sardine_anchovy_sst$sio_sst
	np_sst = sardine_anchovy_sst$np_sst

	# visualize parameters
		par(mfrow=c(2,1))
		par(mar=c(4,4,0.5,0.5))
		plot(sio_sst,type="l",col="mediumseagreen",lwd=2,ylim=range(c(sio_sst,np_sst)),xlab="year",ylab="temperature")
		lines(np_sst,col="lightcoral",lwd=2)

		plot(anchovy,type="l",col="darkorange1",lwd=2,ylim=range(c(anchovy,sardine)),xlab="year",ylab="abundance")
		lines(sardine,col="royalblue1",lwd=2)

	# visualize phase space
		# Z = anchovy
		# plot3d(Z[3:78],Z[2:77],Z[1:76],xlab = "Z", ylab = "Z[i-1]", zlab = "Z[i-2]", type = "l", col="mediumseagreen",lwd=2)

		# Z = np_sst
		# plot3d(Z[3:78],Z[2:77],Z[1:76],xlab = "Z", ylab = "Z[i-1]", zlab = "Z[i-2]", type = "l", col="mediumseagreen",lwd=2)

# EDM analysis
	# package time series
		ts.df = cbind(anchovy,sardine,sio_sst,np_sst)
		colnames(ts.df) = c("anchovy","sardine","sio_sst","np_sst")
		ts.df = as.data.frame(ts.df)

	# calculate dimensionality (simplex)
		max.dim = 5
		simp.loo = apply(ts.df,
						2,
						simplex,
						E=1:max.dim)
		simp.50_50 = apply(ts.df,
						2,
						simplex,
						E=1:max.dim,
						lib=c(1,round(0.5*nrow(ts.df))),
						pred=c(round(0.5*nrow(ts.df))+1,nrow(ts.df)))

		# visualize
			# ts.cols = colorRampPalette(c("gray80","#7d82d3","#ce4c71","#f91741"), space = c("Lab"))(length(simp.loo))
			ts.cols = c("darkorange1","royalblue1","mediumseagreen","lightcoral")

			par(mfrow = c(2,1),mar=c(4,4,0.5,0.5))
			plot(1,1,type="n",xlim=c(1,max.dim),ylim=c(min(sapply(simp.loo,function(x)min(x$rho))),1),xlab="E",ylab="rho")
			for(i in 1:length(simp.loo))
			{
				lines(simp.loo[[i]]$E,simp.loo[[i]]$rho,col=ts.cols[i],lwd=2)
			}
			legend("topright",legend=names(simp.loo),lwd=2,col=ts.cols)

			plot(1,1,type="n",xlim=c(1,max.dim),ylim=c(min(sapply(simp.50_50,function(x)min(x$rho))),1),xlab="E",ylab="rho")
			for(i in 1:length(simp.50_50))
			{
				lines(simp.50_50[[i]]$E,simp.50_50[[i]]$rho,col=ts.cols[i],lwd=2)
			}
			legend("topright",legend=names(simp.50_50),lwd=2,col=ts.cols)

		# optimal embedding dimension
			E.loo = sapply(simp.loo,function(x) x$E[which(x$rho == max(x$rho))])
			rho.loo = sapply(simp.loo,function(x) max(x$rho))
			E.50_50 = sapply(simp.50_50,function(x) x$E[which(x$rho == max(x$rho))])
			rho.50_50 = sapply(simp.50_50,function(x) max(x$rho))
		
	# determine non-linearity (s-maps)
		theta.vec = c(0, 1e-04, 3e-04, 0.001,0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10)
		
		smap.list = as.list(rep(NA,ncol(ts.df)))
		names(smap.list) = colnames(ts.df)

		for(i in 1:ncol(ts.df))
		{
			smap.list[[i]] = s_map(ts.df[,i],E=E.loo[i],theta=theta.vec)
		}

		# visualize
			plot(1,1,type="n",xlim=range(theta.vec),ylim=c(min(sapply(smap.list,function(x)min(x$rho))),1),xlab="theta",ylab="rho")
			for(i in 1:length(smap.list))
			{
				lines(smap.list[[i]]$theta,smap.list[[i]]$rho,col=ts.cols[i],lwd=2)
			}
			legend("topright",legend=names(smap.list),lwd=2,col=ts.cols)

		# optimal theta
			theta.opt = sapply(smap.list,function(x) x$theta[which(x$rho == max(x$rho))])

		# bootstrap for significance test of non-linearity
			# nlboot.list = as.list(rep(NA,ncol(ts.df)))
			# names(nlboot.list) = colnames(ts.df)

			# for(i in 1:ncol(ts.df))
			# {
			# 	nlboot.list[[i]] = test_nonlinearity(ts.df[,i],E=E.loo[i],theta=theta.vec,method="random_shuffle",num_surr=50)
			# }

			# nl.pval = sapply(nlboot.list,function(x)x$delta_mae_p_value)

	# convergent cross mapping (CCM)
	# apply CCM across all possible combinations
		ts.combos = expand.grid(index_1 = 1:ncol(ts.df),index_2 = 1:ncol(ts.df))
		# remove self CCM combinations (self predictability assessed as max rho from optimal embedding)
			ts.combos = ts.combos[-which(ts.combos$index_1 == ts.combos$index_2),]
		
		CCM.wrapper = function(x)
		# pass one row of ts.combos
		{
			# identify ts
				ts.1 = ts.df[,x[[1]]]
				ts.2 = ts.df[,x[[2]]]

				# concatenate
					DF = cbind(ts.1, ts.2)
					colnames(DF) = c("ts.1","ts.2")
					n.ts = nrow(DF)
					

				# identify optimal E use it for CCM
					optim.E.df = as.data.frame(matrix(NA,nrow=max.dim,ncol=11))
					colnames(optim.E.df) = c("E","tau","tp","num_neighbors","lib_column","target_column","lib_size","num_pred","rho","mae","rmse")
					for(e in 1:max.dim)
					{
						optim.E.df[e,] = ccm_means(ccm(DF, lib_column = "ts.1", target_column = "ts.2", E = e, lib_sizes = n.ts, tp = -1,RNGseed=777,silent=TRUE))
					}
					
					optim.E = optim.E.df$E[which.max(optim.E.df$rho)]

				if(length(optim.E) == 0)
				{
					return(c(NA,NA))
				} else {
						# co-predict (ts.1 to ts.2)
							ccm.pair = ccm_means(ccm(DF, lib_column = "ts.1", target_column = "ts.2", E = optim.E, lib_sizes = n.ts, tp = 0,RNGseed=777,silent=TRUE))
							ccm.rho = ccm.pair$rho

							return(c(optim.E,ccm.rho))
				}
		}

	# parallel
		# install.packages("snow")
		# install.packages("snowfall")
		library(snow)
		library(snowfall)
	cores =  4  #  Change as required

	A = proc.time()
		sfInit(parallel=TRUE, cpus=cores, type="SOCK")
		sfExport(list=c("ts.df","max.dim"))
		sfLibrary(rEDM)
		sfClusterSetupRNG(seed=777)
		CCM.output = t(sfApply(as.matrix(ts.combos),1,CCM.wrapper))
		sfRemoveAll()
		sfStop()
	B = proc.time()
	B-A

	CCM.output.DF = data.frame(ts_1 = colnames(ts.df)[ts.combos$index_1],
	ts_2 = colnames(ts.df)[ts.combos$index_2],
	E.dim = CCM.output[,1],
	rho = CCM.output[,2])

	# take a closer look
	ts1 = "anchovy"
	ts2 = "np_sst"					
	optim.E = 2

	anchovy.xmap.np_sst = ccm_means(ccm(ts.df, lib_column = ts1, target_column = ts2, E = optim.E, lib_sizes = round(seq(from=1,to=nrow(ts.df),length.out=20)), tp = 0,RNGseed=777,silent=TRUE))

	ts1 = "np_sst"
	ts2 = "anchovy"
	optim.E = 3

	np_sst.xmap.anchovy = ccm_means(ccm(ts.df, lib_column = ts1, target_column = ts2, E = optim.E, lib_sizes = round(seq(from=1,to=nrow(ts.df),length.out=20)), tp = 0,RNGseed=777,silent=TRUE))

	plot(anchovy.xmap.np_sst$lib_size, pmax(0, anchovy.xmap.np_sst$rho), type = "l", col = "red",lwd=2, xlab = "Library Size", 
    ylab = "Cross Map Skill (rho)", ylim = c(0, 0.3))
    lines(np_sst.xmap.anchovy$lib_size, pmax(0, np_sst.xmap.anchovy$rho),col="blue",lwd=2)
    legend("topright",legend=c("anchovy.xmap.np_sst","np_sst.xmap.anchovy"),lwd=2,col=c("red","blue"))
