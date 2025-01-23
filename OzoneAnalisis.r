library(plot3D)

source("JDR.BT.r")

datosO3 <- read.table("data-O3.txt",header=TRUE)
names(datosO3)
attach(datosO3)

dim(datosO3)
d03=datosO3[20060300<datosO3$date & datosO3$date<20060400,] #March,2006
names(d03)
data1=d03[,4:10]
dim(data1)
names(data1)

install.packages("aTSA")
library(aTSA)

adf.test(data1[,1])
adf.test(data1[,2])
adf.test(data1[,3])
adf.test(data1[,4])
adf.test(data1[,5])
adf.test(data1[,6])
adf.test(data1[,7])

k=dim(data1)[2] # Number of Stations

J <- matrix(NA,k,k)
diag(J) = 0

i=6; j=7  # Raro que no funcione el for
J[i,j] = JDR.BT(data1[,i], data1[,j])

J <- round(J,3)
colnames(J) <- rownames(J) <- names(datosO3)[-c(1:3)]
J