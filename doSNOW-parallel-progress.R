#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#			EXAMPLE CODE OF FOREACH PARALLEL PROCESSING
#				-doSNOW style w/ progress bar
#		
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


###-----------------------------------------------------
#		Packages
###-----------------------------------------------------
library(doSNOW)
library(foreach)

###-----------------------------------------------------
#		Data
###-----------------------------------------------------
mydata <- matrix(rnorm(8000*500), ncol=500)
mydata[sample.int(8000*500, size=4000)] <- NA

###-----------------------------------------------------
#		Register Cores
###-----------------------------------------------------
cores <- parallel::detectCores()-1
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)

###-----------------------------------------------------
#		Make progress bar
###-----------------------------------------------------
pb <- txtProgressBar(min=1, max=8000, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

###-----------------------------------------------------
#		FOREACH LOOP
###-----------------------------------------------------
st <- Sys.time()
result <- foreach(i=1:8000, .options.snow=opts, .combine='rbind') %dopar% {
		v <- which(is.na(mydata[i,]))
		if(length(v)>0){	
			for(j in 1:length(v)){
				if(v[j] != 1 & v[j] != ncol(mydata)){
					mydata[i,v[j]] <- mean(mydata[i,c(v[j]-1,v[j]+1)])
				}else if(v[j] == 1){
					mydata[i,v[j]] <- mydata[i,v[j]+1]
				}else if(v[j] == ncol(mydata)){
					mydata[i,v[j]] <- mydata[i,v[j]-1]
				}
			}
		}
		mydata[i,]
}
close(pb)
stopCluster(cl)
end <- Sys.time()-st
print(end)




