library(LongMemoryTS)
library(entropy)
library(ggplot2)
library(plotrix)
library(GA)

source("JDR.BT.r")
library(plot3D)

# A1 vs A2

n = 500 # 100, 500, 1000
A <- seq(1,10,0.1)
J <- matrix(NA,length(A),length(A))

for(i in 1:length(A)) for(j in 1:length(A)){

	X <- matrix(NA, 2, n)

	for(k in 1:n){
		X[1,k] = A[i]*sin(2*pi*k + 3*pi)
		X[2,k] = A[j]*sin(2*pi*k + 3*pi)		
	}

	J[i,j] = JDR.BT(X[1,], X[2,])
}

#round(J, 4)

image2D(z=J*1e+25, y=A, x=A, main=expression(phi~"="~3*pi), 
clab=expression(10^25*JDR[alpha](x,y)~"="),
xlab=expression(A[1]), ylab=expression(A[2]))  #

# phi1 vs phi2

n = 500 # 100, 500, 1000
A = 10
r <- seq(0,3,0.025)
phi <- r*pi
J <- matrix(NA,length(phi),length(phi))

for(i in 1:length(phi)) for(j in 1:length(phi)){

	X <- matrix(NA, 2, n)

	for(k in 1:n){
		X[1,k] = A*sin(2*pi*k + phi[i])
		X[2,k] = A*sin(2*pi*k + phi[j])		
	}

	J[i,j] = JDR.BT(X[1,], X[2,])
}

#round(J, 4)

image2D(z=J, y=phi, x=phi, main=expression(A~"="~10), 
clab=expression(JDR[alpha](x,y)~"="),
xlab=expression(phi[1]), ylab=expression(phi[2]))


