
VDR <- function(d1,d2,sigma1,sigma2,sigma12) {
	d12 = 0.5*(d1+d2)
	g = sigma12*exp(Re(1i*pi*(d2-d1)/2)) / (2*pi)
	faux <- function(lambda) log(sigma1*lambda^(-2*d1) + sigma2*lambda^(-2*d2) - 2*(2*pi)*g*lambda^(-2*d12))
	R = integrate(faux, 0, pi)$value / pi
	return(R)
}

## FN process

#install.packages("SimDesign")
library(LongMemoryTS)
library(SimDesign)
library(plot3D)
source("VDR.Est.r")

fft.aux <- function(x){
	n = length(x)
	k <- c(1:n)
	f <- rep(NA,n)
	
	for(h in 1:n) {
		f[h] = sum(  x[k]*exp(-1i*2*pi*(k-1)*(h-1)/n)  )
	}

	return(Re(f))
}

fdiff <- function (x, d){  # Jensen-Nielsen method
    iT <- length(x)
    np2 <- nextn(2 * iT - 1, 2)
    k <- 1:(iT - 1)
    b <- c(1, cumprod((k - d - 1)/k))
    dx <- fft(fft(c(b, rep(0, np2 - iT))) * fft(c(x, rep(0, np2 - iT))), 
                  inverse = T)/np2
    return(Re(dx[1:iT]))
}

rho = 0.35
Sigma <- diag(2)
Sigma[1,2] = Sigma[2,1] = rho
n = 100  # 50, 100, 200, 400
phi = THETA = 0
d1 <- d2 <- seq(-0.99,0.49,0.01)
b = 250
noise <- rmvnorm((n+b),mean=rep(0,2),sigma=Sigma)
V <- matrix(NA,length(d1),length(d2))

for(i in 1:length(d1)) for(j in 1:length(d2)) {
	D = c(d1[i],d2[j])
	X <- matrix(NA,(n+b),2)
	for(k in 1:2) X[,k] <- fdiff(noise[,k],-D[k])
	V[i,j] = VDR.Est(x=X[(b+1):(b+n),1],y=X[(b+1):(b+n),2])
}

image2D(z=V, y=d2, x=d1, main=expression(n~"= 400"), 
ylab=expression(d[y]), xlab=expression(d[x]))
