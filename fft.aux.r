fft.aux <- function(x){
	n = length(x)
	k <- c(1:n)
	f <- rep(NA,n)
	
	for(h in 1:n) {
		f[h] = sum(  x[k]*exp(-1i*2*pi*(k-1)*(h-1)/n)  )
	}

	return(Re(f))
}

JDR <- function(d1,d2,sigma1,sigma2,sigma12,alpha) {
	alpha.aux = alpha*(1-alpha)
	zeta.d <- function(d) (1-2*d)*(2*pi)^(1+2*d)
	d12 = 0.5*(d1+d2)
	R = alpha.aux*(sigma1^2/zeta.d(d1) + sigma2^2/zeta.d(d2) - 2*sigma12*exp(Re(1i*pi*(d2-d1)/2))/zeta.d(d12)) 
	return(R)
}

data(lynx)
x <- scale(lynx)
length(x)
plot(fft.aux(x), type="l")

lagEstimate <- function(x,k,N=length(x)) (1/N)*sum(x[1:(N-k)]*x[(k+1):N])
Lags <- function(x,kMax) sapply(0:kMax,lagEstimate,x=x)
acsWindowed <- function(x,kMax,Nzero=0){
  rHalf <- c(Lags(x,kMax),rep(0,Nzero))
  c(rev(rHalf[2:length(rHalf)]),rHalf)
}
Textbook2R <- function(x,N=length(x),foldN=ceiling(N/2)) {
  c(x[foldN:N],x[1:(foldN-1)])
}
BartlettWindow <- function(N,n=seq(0,(N-1))) 1-abs((n-(N-1)/2)/((N-1)/2))

rHat    <- Textbook2R(acsWindowed(x,kMax=65))
W      <- Textbook2R(BartlettWindow(length(rHat)))
FFT_W <- fft.aux(rHat*W)
plot(FFT_W, type="l")
abline(h=0,lty=2)

JDR.aux <- function(x){
	n = length(x)
	h <- c(1:n)
	psi = sin(pi*k)/(pi*k)
 	J = x[1] + 2*sum(x[-1]*psi)
	return(J)
}

rHat <- Textbook2R(acsWindowed(x,kMax=65))
Wb <- Textbook2R(BartlettWindow(length(rHat)))
JDR.aux(rHat*Wb)

lagEstimate.xy <- function(x,y,k,N=length(x)) (1/N)*sum(x[1:(N-k)]*y[(k+1):N])
Lags.xy <- function(x,y,kMax) sapply(0:kMax,lagEstimate.xy,x=x,y=y)
acsWindowed.xy <- function(x,y,kMax,Nzero=0){
  rHalf <- c(Lags.xy(x,y,kMax),rep(0,Nzero))
  c(rev(rHalf[2:length(rHalf)]),rHalf)
}

JDR <- function(x,y,alpha=0.5,kMax=length(x)/5){
	beta = alpha*(1-alpha)
	rHat.x    <- Textbook2R(acsWindowed(x,kMax))
	Wb.x      <- Textbook2R(BartlettWindow(length(rHat.x)))
	rHat.y    <- Textbook2R(acsWindowed(y,kMax))
	Wb.y      <- Textbook2R(BartlettWindow(length(rHat.y)))
	rHat.xy    <- Textbook2R(acsWindowed.xy(x,y,kMax))
	Wb.xy      <- Textbook2R(BartlettWindow(length(rHat.xy)))

	J = JDR.aux(rHat.x*Wb.x) + JDR.aux(rHat.y*Wb.y) - 2*JDR.aux(rHat.xy*Wb.xy)
	return(beta*J)
}

## FN process

#install.packages("SimDesign")
library(LongMemoryTS)
library(SimDesign)
library(plot3D)
source("JDR.BT.r")

fdiff <- function (x, d){  # Jensen-Nielsen method
    iT <- length(x)
    np2 <- nextn(2 * iT - 1, 2)
    k <- 1:(iT - 1)
    b <- c(1, cumprod((k - d - 1)/k))
    dx <- fft(fft(c(b, rep(0, np2 - iT))) * fft(c(x, rep(0, np2 - iT))), 
                  inverse = T)/np2
    return(Re(dx[1:iT]))
}

rho = 0.5 # 0, 0.25, 0.5, 1
Sigma <- diag(2)
Sigma[1,2] = Sigma[2,1] = rho
n = 100  # 50, 100, 200, 400
phi = THETA = 0
d1 <- d2 <- seq(-0.99,0.49,0.01)
b = 250
noise <- rmvnorm((n+b),mean=rep(0,2),sigma=Sigma)
J <- matrix(NA,length(d1),length(d2))

for(i in 1:length(d1)) for(j in 1:length(d2)) {
	D = c(d1[i],d2[j])
	X <- matrix(NA,(n+b),2)
	for(k in 1:2) X[,k] <- fdiff(noise[,k],-D[k])
	J[i,j] = JDR.BT(x=X[(b+1):(b+n),1],y=X[(b+1):(b+n),2])
}

image2D(z=J, y=d2, x=d1, main=expression(n~"= 100"), 
ylab=expression(d[y]), xlab=expression(d[x]),
clab=expression(JDR[alpha](x,y)))

### Comparisson of estimation methods

rho = 0.5 # 0, 0.25, 0.5, 1
Sigma <- 1*diag(2)
Sigma[1,2] = Sigma[2,1] = rho
n = 1000  # 50, 100, 200, 400
phi = THETA = 0
d1 = 0.49 # -0.9, -0.5, -0.25, 0, 0.25, 0.49 
d2 <- seq(-0.99,0.49,0.01)
b = 250
noise <- rmvnorm((n+b),mean=rep(0,2),sigma=Sigma)
J <- matrix(NA,2,length(d2))

for(j in 1:length(d2)) {
	D = c(d1,d2[j])
	X <- matrix(NA,(n+b),2)
	for(k in 1:2) X[,k] <- fdiff(noise[,k],-D[k])
	J[1,j] = JDR.BT(x=X[(b+1):(b+n),1],y=X[(b+1):(b+n),2], window="Bartlett")
	J[2,j] = JDR(d1,d2[j],sigma1=Sigma[1,1],sigma2=Sigma[2,2],sigma12=rho,alpha=0.5) # Exact
}

plot(d2, J[1,], type="l", lwd=2, col="red", ylim=c(min(J),max(J)), lty=2, 
ylab="JDR", xlab=expression(d[y]), main=expression(d[x]~"= 0.49"))
lines(d2, J[2,], lwd=2, col="black")
legend("top",c("Exact","BT est."),col=c("black","red"),
lwd=2,lty=c(1,2),bty="n")
