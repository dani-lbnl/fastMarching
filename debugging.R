setwd("/Users/dani/Dropbox/prog/Apps_IDEAL/fastMarching/Rfastmarching/Rfastmarching_july7")
graphics.off();
source('meshgrid.R')
require(fields)
require(PolynomF)

eps = 1.110223e-16; #NO F=0!!!!!
ex = 1

ptm <- proc.time()

#Example: simple
if (ex==1){
Seeds=matrix(c(40,1),2,1,byrow=F)
F = matrix(1,50,50);
source('msfm2d3.R'); t=msfm2d3(F,Seeds); image(t)
windows()
source('msfm2d4.R'); t=msfm2d4(F,Seeds);image(t)
}

if(ex==2){
#F = matrix(1,101,101); #speed function
Seeds = matrix(c(25,25),2,1,byrow=F); #initial/source points 
M = meshgrid(1:51,1:51)
F = sqrt( (M$x-Seeds[1])^2 + (M$y-Seeds[2])^2 )
halfx=trunc(dim(F)[1]/2)
halfy=trunc(dim(F)[2]/2)
F[halfx:dim(F)[1],halfy:dim(F)[2]]=2;
source('msfm2d3.R'); t=msfm2d3(F,Seeds); image(t)
windows()
source('msfm2d4.R'); t=msfm2d4(F,Seeds);image(t)
}

if(ex==5){
#Supersimple
Seeds = matrix(c(3,3),2,1,byrow=F); #initial/source points 
F = matrix(1,5,5);
source('msfm2d3.R'); t=msfm2d3(F,Seeds); image(t)
windows()
source('msfm2d4.R'); t=msfm2d4(F,Seeds);image(t)
}

if(ex==6){
print('Random seeds, constant F') #image on the R abstract
Seeds=matrix(trunc(runif(20,1,50)),nrow=2,ncol=10,byrow=F)
F = matrix(1,50,50);
source('msfm2d3.R'); t=msfm2d3(F,Seeds); image(t)
windows()
source('msfm2d4.R'); t=msfm2d4(F,Seeds);image(t)
}

if(ex==7){
P = matrix(c(51,51),2,1,byrow=F); #for the speed function 
M = meshgrid(1:101,1:101)
F = sqrt( (M$x-P[1])^2 + (M$y-P[2])^2 ) + 1

#Seeds=matrix(c(45:57,45:57),nrow=2,byrow=T)
Seeds=matrix(c(45:57,seq(57,45,by=-1)),nrow=2,byrow=T);

#source('msfm2d3.R'); t=msfm2d3(F,Seeds); image(t)
#rm(t)
#windows()
source('msfm2d4.R'); t=msfm2d4(F,Seeds);
SeedsM=matrix(rep(0,length(F)),dim(F)[1],dim(F)[2]);
for (i in 1:(dim(Seeds)[2]))
{ SeedsM[Seeds[1,i],Seeds[2,i]]=1
  #F[Seeds[1,i],Seeds[2,i]]=1
}
png(filename=paste('example.png',sep=''), width = 800, height =800, pointsize=12 , bg="white")
par(mfrow=c(2,2))
image(SeedsM,main="Seeds",col=tim.colors(20))
image(F,main="Speed function",col=tim.colors(20))
image(t,main="Result of propagation",col=tim.colors(20))
contour(t,main="Isocontours of propagation")
dev.off()
}

if(ex==8){
mmax=101;
Seeds=matrix(trunc(runif(20,2,(mmax-2))),nrow=2,ncol=10,byrow=F)
M = meshgrid(1:mmax,1:mmax)
F = sqrt( (M$x-Seeds[1])^2 + (M$y-Seeds[2])^2 )+ 
source('msfm2d3.R'); t=msfm2d3(F,Seeds); image(t,col=tim.colors(20))
rm(t)
windows()
source('msfm2d4.R'); t=msfm2d4(F,Seeds);image(t)
}


if(ex==9){
nmax=101;
Seeds=matrix(c(51,51),nrow=2,ncol=3,byrow=F)
x=1:nmax;
#fx=exp(-(x-50)^2/16)
fx=c(rep(0,30),rep(1,41),rep(0,30));
F = matrix(rep(fx,nmax),nrow=nmax,ncol=nmax,byrow=F);
F[1:2,]=0

source('msfm2d3.R'); t1=msfm2d3(F,Seeds); 
source('msfm2d4.R'); t2=msfm2d4(F,Seeds);
image(t1,main='Standard',col=tim.colors(10)); #contour(t,add=T)
windows()
image(t2,main='Timer included',,col=tim.colors(10)); #contour(t,add=T)
}


print(proc.time() - ptm)
png(filename=paste('speed.png',sep=''), width = 800, height =800, pointsize=12 , bg="white")
persp(F, theta = 30, phi = 20, col = "lightblue",
           ltheta = 120, shade = 0.75, ticktype = "detailed",
           xlab = "X", ylab = "Y", zlab = "Moving ability")
dev.off()
