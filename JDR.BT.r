JDR.BT <- function(x, y, alpha=0.5, kMax=length(x)/5, window="Bartlett"){
	
	beta = alpha*(1-alpha)

	# Auxiliary functions

	lagEstimate <- function(x,k,N=length(x)) (1/N)*sum(x[1:(N-k)]*x[(k+1):N], na.rm=TRUE)

	Lags <- function(x,kMax) sapply(0:kMax,lagEstimate,x=x)

	acsWindowed <- function(x,kMax,Nzero=0){
  			rHalf <- c(Lags(x,kMax),rep(0,Nzero))
  			c(rev(rHalf[2:length(rHalf)]),rHalf)
	}

	Textbook2R <- function(x,N=length(x),foldN=ceiling(N/2)) c(x[foldN:N],x[1:(foldN-1)])

	RecWindow <- function(N,n=seq(0,(N-1))) rep(1,N)
	BartlettWindow <- function(N,n=seq(0,(N-1))) 1-abs((n-(N-1)/2)/((N-1)/2))

	lagEstimate.xy <- function(x,y,k,N=length(x)) (1/N)*sum(x[1:(N-k)]*y[(k+1):N], na.rm=TRUE)

	Lags.xy <- function(x,y,kMax) sapply(0:kMax,lagEstimate.xy,x=x,y=y)

	acsWindowed.xy <- function(x,y,kMax,Nzero=0){
  			rHalf <- c(Lags.xy(x,y,kMax),rep(0,Nzero))
  			c(rev(rHalf[2:length(rHalf)]),rHalf)
	}

	JDR.aux <- function(x){
		n = length(x); h <- c(1:(n-1))
		psi = sin(pi*h)/(pi*h)
 		J = x[1] + 2*sum(x[-1]*psi, na.rm=TRUE)
		return(J)
	}

	rHat.x <- Textbook2R(acsWindowed(x,kMax))
	if(window=="Rectangular") Wb.x <- Textbook2R(RecWindow(length(rHat.x)))
	if(window=="Bartlett") Wb.x <- Textbook2R(BartlettWindow(length(rHat.x)))

	rHat.y <- Textbook2R(acsWindowed(y,kMax))
	if(window=="Rectangular") Wb.y <- Textbook2R(RecWindow(length(rHat.y)))
	if(window=="Bartlett") Wb.y <- Textbook2R(BartlettWindow(length(rHat.y)))

	rHat.xy <- Textbook2R(acsWindowed.xy(x,y,kMax))
	if(window=="Rectangular") Wb.xy <- Textbook2R(RecWindow(length(rHat.xy)))
	if(window=="Bartlett") Wb.xy <- Textbook2R(BartlettWindow(length(rHat.xy)))

	J = JDR.aux(rHat.x*Wb.x) + JDR.aux(rHat.y*Wb.y) - 2*JDR.aux(rHat.xy*Wb.xy)
	return(beta*J)
}