source("use_for_all/matrix_generator.R")
seed <- 1:400

#Dimension
p <- 200 #400
load(file=paste("L0_Gaussian_band_p",p,"_supportset_of_normalplot.RData",sep=""))
#Interset of support set of n=200,400,800
ground_support <- which(is.na(omega_support))

#Sample size
n <- 200 #400,800

#Z score 
z_100 <- {}
for(i in seed){
  dir <- paste("L1_Gaussian_normalplot_output/L1_desparsified_Gaussian_band_n",n,"_p",p,"_seed",i,"_normalplot.RData",sep="")
  load(file=dir)
  omega3 <- matrix(1,nrow=p,ncol=p)
  diag(omega3) <- NA
  omega3_sigma_c <- matrix(omega3[which(!is.na(omega3))],nrow=p-1)
  omega3_sigma_c[Ahat] <- NA
  omega3_support<-matrix_generator3(diag(omega3),omega3_sigma_c)
  
  Ahat1 <- which(is.na(omega3_support))
  Ahat2 <- which(Ahat1 %in% ground_support)
  
  z_100 <- cbind(z_100,z[Ahat2])
  rm(Ahat,z,sdar_F_sigma,NS_F_sigma,NS_T,data)
}

u <- seq(-3,3,by=0.01)

hist(z_100[3,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.45),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(1,2)",breaks=seq(-20,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_100[3,]),col="blue3",lwd=2.5,lty=1)

hist(z_100[5,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.45),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(3,2)",breaks=seq(-20,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_100[5,]),col="blue3",lwd=2.5,lty=1)

hist(z_100[6,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.45),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(2,3)",breaks=seq(-20,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_100[6,]),col="blue3",lwd=2.5,lty=1)

hist(z_100[8,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.45),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(4,3)",breaks=seq(-20,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_100[8,]),col="blue3",lwd=2.5,lty=1)



save.image(paste("L1_desparsified_Gaussian_band_n",n,"_p",p,"_normalplot.RData",sep=""))