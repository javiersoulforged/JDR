
# Entropy rates

SER <- function(d,s2){ # Shannon entropy rate
	Cx = s2*sin(pi*(1-d/2))*gamma(3-d)
	0.5*log(2*pi*exp(1)*Cx) + (d-1)*(log(pi)-1)
}

RER <- function(d,s2){ # quadratic Rényi entropy rate
	Cx = s2*sin(pi*(1-d/2))*gamma(3-d)
	0.5*log(4*pi*Cx) + (d-1)*(log(pi)-1)
}

d <- seq(0.005,2,0.01)

plot(d, SER(d,1), type="l", lwd=2, ylim=c(-1,2.25), ylab="Entropy rates")
lines(d, RER(d,1), type="l", lwd=2, lty=2)
lines(d, SER(d,2), type="l", col="blue", lwd=2)
lines(d, RER(d,2), type="l", col="blue", lwd=2, lty=2)
lines(d, SER(d,3), type="l", col="red", lwd=2)
lines(d, RER(d,3), type="l", col="red", lwd=2, lty=2)
lines(d, SER(d,4), type="l", col="purple", lwd=2)
lines(d, RER(d,4), type="l", col="purple", lwd=2, lty=2)
legend("bottom", c(expression("SER"~"("~sigma^2~"= 1)"), 
expression("SER"~"("~sigma^2~"= 2)"),
expression("SER"~"("~sigma^2~"= 3)"),
expression("SER"~"("~sigma^2~"= 4)"),
expression("RER"~"("~sigma^2~"= 1)"), 
expression("RER"~"("~sigma^2~"= 2)"),
expression("RER"~"("~sigma^2~"= 3)"),
expression("RER"~"("~sigma^2~"= 4)")), 
col=c("black","blue","red","purple","black","blue","red","purple"), 
lty=c(1,1,1,1,2,2,2,2),lwd=2, bty="n", ncol=2)

# MI rates

library(plot3D)

MIR <- function(rho) -0.5*log(1-rho^2)

rho <- seq(-1,1,0.001)

plot(rho, MIR(rho), ylab="MIR", xlab=expression(rho))

# Distance rates

JDR <- function(d1,d2,sigma1,sigma2,sigma12,alpha) {
	alpha.aux = alpha*(1-alpha)
	zeta.d <- function(d) (1-2*d)*(2*pi)^(1+2*d)
	d12 = 0.5*(d1+d2)
	R = alpha.aux*(sigma1^2/zeta.d(d1) + sigma2^2/zeta.d(d2) - 2*sigma12*exp(Re(1i*pi*(d2-d1)/2))/zeta.d(d12)) 
	return(R)
}


JVDR <- function(d1,d2,sigma1,sigma2,sigma12) {
	d12 = 0.5*(d1+d2)
	g = sigma12*exp(Re(1i*pi*(d2-d1)/2)) / (2*pi)
	faux <- function(lambda) log(sigma1*lambda^(-2*d1) + sigma2*lambda^(-2*d2) - 2*(2*pi)*g*lambda^(-2*d12))
	R = integrate(faux, 0, pi)$value / pi
	return(R)
}

d1 <- d2 <- seq(-0.99,0.49,0.01)
J <- matrix(NA,length(d1),length(d2))

for(i in 1:length(d1)) for(j in 1:length(d2)){
	J[i,j] = JDR(d1[i],d2[j],1,1,0.35,0.5)
}

image2D(z=J, y=d2, x=d1, xlab=expression(d[x]), ylab=expression(d[y]))

g <- function(lambda,d1,d2,sigma1,sigma2,sigma12) {
	d12 = 0.5*(d1+d2)
	sigma1*lambda^(-d1) + sigma2*lambda^(-d2) - 2*sigma12*exp(Re(1i*pi*(d2-d1)/2))*lambda^(-2*d12)
}

lambda <- seq(0,pi,0.01)

plot(lambda, g(lambda, 0.25, 0.1, 2, 2, 1))



d1 <- d2 <- seq(-0.99,0.48,0.01)
J <- matrix(NA,length(d1),length(d2))

for(i in 1:length(d1)) for(j in 1:length(d2)){
	J[i,j] = JVDR(d1[i],d2[j],1,1,0.35)
}

image2D(z=J, y=d2, x=d1, xlab=expression(d[x]), ylab=expression(d[y]))

