source('meshgrid.R')
require(PolynomF)
#Example: simple
Seeds=matrix(trunc(runif(10,1,20)),2,50,byrow=F)
F = matrix(1,20,20);
source('msfm2d4.R'); t=msfm2d4(F,Seeds);image(t)
#source('msfm2d3.R'); t=msfm2d3(F,Seeds);image(t)

#F = matrix(1,101,101); #speed function
Seeds = matrix(c(51,51),2,1,byrow=F); #initial/source points 
M = meshgrid(1:101,1:101)
F = sqrt( (M$x-Seeds[1])^2 + (M$y-Seeds[2])^2 )

halfx=trunc(dim(F)[1]/2)
halfy=trunc(dim(F)[2]/2)
F[halfx:dim(F)[1],halfy:dim(F)[2]]=2;



source('msfm2d3.R'); 
t=msfm2d3(F,Seeds)
image(t)

#Supersimple
Seeds = matrix(c(3,3),2,1,byrow=F); #initial/source points 
F = matrix(1,5,5);
source('msfm2d3.R'); t=msfm2d3(F,Seeds);image(t)

