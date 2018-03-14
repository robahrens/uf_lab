source("read.admb.R")
require(KernSmooth)
"panel.cor" <- function(x, y, digits=2, prefix="", cex.cor,prior,A,B,...)
{
    usr <- par("usr"); on.exit(par(usr))
    v=cor(x[1:length(x)/2], y[1:length(x)/2])
	ymin=usr[3]
	ymax=usr[4]
	yrng=ymax-ymin
	ymin1=usr[3]+.1*yrng
	ymax1=usr[4]-.2*yrng
	#beerclr="darkgoldenrod2"
	#bubbleclr="darkgoldenrod1"
	yrng1=ymax1-ymin1
	ymax2=ymin1+abs(v)*yrng1
	xmid=(usr[2]+usr[1])/2
	ymid=(ymax1+ymin1)/2
	xrng=(usr[2]-usr[1])
	#xpoly=c(xmid-.15*xrng,xmid-(.15+abs(v)*.1)*xrng,xmid+(.15+abs(v)*.1)*xrng,xmid+.15*xrng,xmid-.15*xrng)
	#ypoly=c(ymin1,ymax2,ymax2,ymin1,ymin1)
	#polygon(xpoly,ypoly,col=beerclr,border=beerclr)
	#nbub=round(500*abs(v),0)
	#bubblex=runif(nbub,xmid-(.15+abs(v*.95)*.1)*xrng,xmid+(.15+abs(v*.95)*.1)*xrng)
	#bubbley=runif(nbub,ymax2-.02*yrng1,ymax2+.02*yrng1)
	#points(bubblex,bubbley,pch=21,col = bubbleclr, bg = "white",cex=runif(nbub,1,2))
	#points(c(xmid-.15*xrng,xmid+.15*xrng),c(ymin1,ymin1),type="l",lwd=4)
	#points(c(xmid-.15*xrng,xmid-.25*xrng),c(ymin1,ymax1),type="l",lwd=4)
	#points(c(xmid+.15*xrng,xmid+.25*xrng),c(ymin1,ymax1),type="l",lwd=4)
	if(v<0){
		text(xmid,ymid,labels=c(paste("-",round(abs(v),2),sep="")),cex=(.8+abs(v)*1))
	}else{
		text(xmid,ymid,labels=c(paste("+",round(abs(v),2),sep="")),cex=(.8+abs(v)*1))
	
	}
		
}

"panel.hist" <- function(x,prior,A,B,...)
{
    usr <- par("usr"); on.exit(par(usr))
	print(usr)    
    par(usr = c(usr[1:2], 0, 1.5) )
    ind=par("mfg")[1]
    d=density(A[,ind],adjust=2)
	d$y=d$y/max(d$y)
 	lines(d$x,d$y,lwd=2)
	#rug(sample(x,200),col="black")
	box()
	if(prior[ind]>0)
	{
			dd=density(B[,ind],adjust=2)
			dd$y=dd$y/max(dd$y)
	    	lines(dd$x,dd$y,type="l",lwd=2,col="grey")
	}
}


"panel.contour" <- function(x, y,bw=25,prior,A,B,...)
{
    	usr <- par("usr"); on.exit(par(usr))
	    ind=par("mfg")[2]
	    ind2=par("mfg")[1]
		bw=15
		bwx=(max(A[,ind])-min(A[,ind]))/bw; bwy=(max(A[,ind2])-min(A[,ind2]))/bw
		est <- bkde2D(cbind(A[,ind],A[,ind2]),bandwidth=c(bwx,bwy),gridsize=c(81, 81))
		est$fhat=est$fhat/max(est$fhat)
		#text(max(x),max(y),labels="D",adj=c(1,1))
		lvs=c(0.05,0.5,0.95)
		maxct=max(lvs)
		nlvs=length(lvs)
		#thelines=contourLines(est$x1,est$x2,est$fhat,levels=lvs)
		#polygon(thelines[[nlvs-3]]$x,thelines[[nlvs-3]]$y,col="white",border="black",lwd=1)
		#polygon(thelines[[nlvs-2]]$x,thelines[[nlvs-2]]$y,col="snow",border="black",lwd=2)
		#polygon(thelines[[nlvs-1]]$x,thelines[[nlvs-1]]$y,col="yellow",border="yellow2",lwd=3)
		#polygon(thelines[[nlvs]]$x,thelines[[nlvs]]$y,col="lightyellow",border="yellow",lwd=1)
		#xi=sample(1:length(x),round(.1*length(x),0))
		#points(x[xi],y[xi],pch=".",col=grey(0:10/10))
		contour(est$x1,est$x2,est$fhat,drawlabels=T,add=T,levels=lvs,lty=1,lwd=1,labcex= 0.6)
		#Add salt and pepper
	box()
}

"plot.pairs"=function(A,txt="",label,prior,mome1,mome2)
{
    npar=length(A[1,])
	nrec=length(A[,1])
	B=matrix(0,nrec,npar)
	colnames(B)=colnames(A)
	for(ind in 1:npar)
	{
	if(prior[ind]>0)
	{
	 	if(prior[ind]==1) B[,ind]=rlnorm(nrec,mome1[ind],mome2[ind])
		if(prior[ind]==2) B[,ind]=sqrt(1./rgamma(nrec,mome1[ind],mome2[ind]))
		if(prior[ind]==3) B[,ind]=rbeta(1000,mome1[ind],mome2[ind])
	}else{
		B[,ind]=A[,ind]
	}
	}
	C=rbind(A,B)
	pairs(C,lower.panel=panel.contour,upper.panel=panel.cor,diag.panel=panel.hist,gap=0,labels=label,prior=prior,A=A,B=B)
	mtext(txt,3,outer=T,line=-1,cex=.85)
}
