source("matrix_generator.R")
seed <- 1:400

#Sample size
n <- 200
#Dimension
p <- 200 #400

Ahat_ground_n200 <- c(1:((p-1)*p)) 
for(i in seed){
  dir <- paste("L0_subGaussian_normalplot_output/L0_desparsified_unbias_subGaussian_band_n",n,"_p",p,"_seed",i,"_normalplot.RData",sep="")
  load(file=dir)
  Ahat_ground_n200 <- intersect(Ahat_ground_n200,Ahat)
  rm(Ahat,z_unbias,z_desparsified,sdar_F_sigma,sdar_T,data)
}

n <- 400
Ahat_ground_n400 <- c(1:((p-1)*p)) 
for(i in seed){
  dir <- paste("L0_subGaussian_normalplot_output/L0_desparsified_unbias_subGaussian_band_n",n,"_p",p,"_seed",i,"_normalplot.RData",sep="")
  load(file=dir)
  Ahat_ground_n400 <- intersect(Ahat_ground_n400,Ahat)
  rm(Ahat,z_unbias,z_desparsified,sdar_F_sigma,sdar_T,data)
}

n <- 800
Ahat_ground_n800 <- c(1:((p-1)*p)) 
for(i in seed){
  dir <- paste("L0_subGaussian_normalplot_output/L0_desparsified_unbias_subGaussian_band_n",n,"_p",p,"_seed",i,"_normalplot.RData",sep="")
  load(file=dir)
  Ahat_ground_n800 <- intersect(Ahat_ground_n800,Ahat)
  rm(Ahat,z_unbias,z_desparsified,sdar_F_sigma,sdar_T,data)
}

#Interset of support set of n=200,400,800
Ahat_ground <- intersect(intersect(Ahat_ground_n200,Ahat_ground_n400),Ahat_ground_n800)
omega <- matrix(1,nrow=p,ncol=p)
diag(omega) <- NA
omega_sigma_c <- matrix(omega[which(!is.na(omega))],nrow=p-1)
omega_sigma_c[Ahat_ground] <- NA
omega_support<-matrix_generator3(diag(omega),omega_sigma_c)

save(omega_support,file=paste("L0_subGaussian_band_p",p,"_supportset_of_normalplot.RData",sep=""))


p <- 200 #400
load(file=paste("L0_subGaussian_band_p",p,"_supportset_of_normalplot.RData",sep=""))
ground_support <- which(is.na(omega_support))
source("matrix_generator.R")
seed <- 1:400

n <- 200 #400,800

#Z score of desparsified estimator
z_desparsified100 <- {}
#Z score of unbias estimator
z_unbias100 <- {}
for(i in seed){
  dir <- paste("L0_subGaussian_normalplot_output/L0_desparsified_unbias_subGaussian_band_n",n,"_p",p,"_seed",i,"_normalplot.RData",sep="")
  load(file=dir)
  omega3 <- matrix(1,nrow=p,ncol=p)
  diag(omega3) <- NA
  omega3_sigma_c <- matrix(omega3[which(!is.na(omega3))],nrow=p-1)
  omega3_sigma_c[Ahat] <- NA
  omega3_support<-matrix_generator3(diag(omega3),omega3_sigma_c)
  
  Ahat1 <- which(is.na(omega3_support))
  #For every replication, find index of support set in intersect of support set of Omega_hat_US under 100 replications
  Ahat2 <- which(Ahat1 %in% ground_support)
  
  z_desparsified100 <- cbind(z_desparsified100,z_desparsified[Ahat2])
  z_unbias100 <- cbind(z_unbias100,z_unbias[Ahat2])
  rm(Ahat,z_unbias,z_desparsified,sdar_F_sigma,sdar_T,data)
}

#Draw normal plot
u <- seq(-3,3,by=0.01)

hist(z_unbias100[3,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.48),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(1,2)",breaks=seq(-10,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_unbias100[3,]),col="blue3",lwd=2.5,lty=1)

hist(z_unbias100[5,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.48),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(3,2)",breaks=seq(-10,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_unbias100[5,]),col="blue3",lwd=2.5,lty=1)

hist(z_unbias100[6,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.48),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(2,3)",breaks=seq(-10,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_unbias100[6,]),col="blue3",lwd=2.5,lty=1)

hist(z_unbias100[8,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.48),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(4,3)",breaks=seq(-10,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_unbias100[8,]),col="blue3",lwd=2.5,lty=1)



hist(z_desparsified100[3,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.48),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(1,2)",breaks=seq(-10,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_desparsified100[3,]),col="blue3",lwd=2.5,lty=1)

hist(z_desparsified100[5,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.48),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(3,2)",breaks=seq(-10,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_desparsified100[5,]),col="blue3",lwd=2.5,lty=1)

hist(z_desparsified100[6,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.48),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(2,3)",breaks=seq(-10,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_desparsified100[6,]),col="blue3",lwd=2.5,lty=1)

hist(z_desparsified100[8,],freq=FALSE,xlim=c(-3,3),ylim=c(0,0.48),
     col="cornflowerblue",main="",ylab="",xlab="(i,j)=(4,3)",breaks=seq(-10,10,0.5),cex.lab=1.8,cex.axis=1.5,yaxt="n")
lines(u,dnorm(u),lwd=2,col="red")
abline(v=0,col="red",lwd=2.5,lty=1)
abline(v=mean(z_desparsified100[8,]),col="blue3",lwd=2.5,lty=1)

save.image(paste("L0_desparsified_unbias_subGaussian_band_n",n,"_p",p,"_normalplot.RData",sep=""))
