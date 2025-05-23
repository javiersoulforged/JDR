VDR.Est <- function(x, y, tau=length(x)/5, Bartlett=TRUE){
	z <- x - y
	lags <- 0:tau
	gamma.z <- acf(z, lag.max=tau, type="covariance", plot = FALSE)$acf		
	rho.h <- acf(z, lag.max=tau, type="correlation", plot = FALSE)$acf

	sinc <- function(h) sin(pi*h)/(pi*h)
	window.bar <- function(h,tau=tau) 1 - abs(h)/tau
	
	if(Bartlett==FALSE) V.aux =  sum(rho.h[-1]*sinc(lags[-1]))
	if(Bartlett==TRUE) V.aux = sum(window.bar(lags[-1],tau)*rho.h[-1]*sinc(lags[-1]))

	V = log(gamma.z[1]) + V.aux

	return(V)
}