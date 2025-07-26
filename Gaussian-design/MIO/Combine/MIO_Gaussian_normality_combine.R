library(huge)
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
seed <- scan("seed.txt",what=integer())
n <- 200 #400
p <- 200
Ahat_ground <- c(1:((p-1)*p)) 
z_desparsified100 <- {}
l_desparsified100 <- {}
ltrue_desparsified100 <- {}
cov_desparsified100 <- {}
for(i in seed){
  dir <- paste("MIO_desparsified_unbias_Gaussian_normality/MIO_desparsified_unbias_Gaussian_band_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
  load(file=dir)
  Ahat_ground <- intersect(Ahat_ground,Ahat)
  z_desparsified100 <- cbind(z_desparsified100,z_desparsified)
  l_desparsified100 <- cbind(l_desparsified100,l_desparsified)
  ltrue_desparsified100 <- cbind(ltrue_desparsified100,ltrue_desparsified)
  cov_desparsified100 <- cbind(cov_desparsified100,cov_desparsified)
  rm(Ahat,z_desparsified,l_desparsified,ltrue_desparsified,cov_desparsified)
}
omega <- matrix(1,nrow=p,ncol=p)
diag(omega) <- NA
omega_sigma_c <- matrix(omega[which(!is.na(omega))],nrow=p-1)
omega_sigma_c[Ahat_ground] <- NA
omega_support<-matrix_generator3(diag(omega),omega_sigma_c)

#Generate Band graph
set.seed(-1)
omega2 <- omega_generator1(p,1,0.5,0.3)
Sigma2 <- solve(omega2)

#Generate Random graph
# graph= "random"
# set.seed(-1)
# generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, prob=4/p)
# omega2 <- generator$omega
# Sigma2 <- generator$sigma

#Generate Hub graph
# graph= "hub"
# set.seed(-1)
# generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
# omega2 <- generator$omega
# Sigma2 <- generator$sigma

#Generate Cluster graph
# graph= "cluster"
# set.seed(-1)
# generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
# omega2 <- generator$omega
# Sigma2 <- generator$sigma

#True support set : S_Omega
S0 <- which(abs(omega2)>0.00001)
#Complementary set of true support set : S_Omega_C
S0C <- which(abs(omega2)<=0.00001)
#Intersect of support set of Omega_hat_US under 100 replications : S_hat_Omega
ground_support <- which(is.na(omega_support))

#For every (i,j), compute mean of asymptotic normality results over 100 replications
z_S0_ij_mean <- abs(apply(z_desparsified100[S0,],1,mean))
z_S0C_ij_mean <- abs(apply(z_desparsified100[S0C,],1,mean))
z_ground_ij_mean <- abs(apply(z_desparsified100[ground_support,],1,mean))
z_ij_mean <- abs(apply(z_desparsified100,1,mean))
z_S0_ij_sd <- apply(z_desparsified100[S0,],1,sd)
z_S0C_ij_sd <- apply(z_desparsified100[S0C,],1,sd)
z_ground_ij_sd <- apply(z_desparsified100[ground_support,],1,sd)
z_ij_sd <- apply(z_desparsified100,1,sd)
l_S0_ij_mean <- apply(l_desparsified100[S0,],1,mean)
l_S0C_ij_mean <- apply(l_desparsified100[S0C,],1,mean)
l_ground_ij_mean <- apply(l_desparsified100[ground_support,],1,mean)
l_ij_mean <- apply(l_desparsified100,1,mean)
ltrue_S0_ij_mean <- apply(ltrue_desparsified100[S0,],1,mean)
ltrue_S0C_ij_mean <- apply(ltrue_desparsified100[S0C,],1,mean)
ltrue_ground_ij_mean <- apply(ltrue_desparsified100[ground_support,],1,mean)
ltrue_ij_mean <- apply(ltrue_desparsified100,1,mean)
cov_S0_ij_mean <- apply(cov_desparsified100[S0,],1,mean)
cov_S0C_ij_mean <- apply(cov_desparsified100[S0C,],1,mean)
cov_ground_ij_mean <- apply(cov_desparsified100[ground_support,],1,mean)
cov_ij_mean <- apply(cov_desparsified100,1,mean)

#compute mean of asymptotic normality results over S_Omega, S_Omega_C and S_hat_Omega
coverage_rate_S0_mean <- mean(cov_S0_ij_mean)
coverage_rate_S0C_mean <- mean(cov_S0C_ij_mean)
coverage_rate_ground_mean <- mean(cov_ground_ij_mean)
coverage_rate_mean <- mean(cov_ij_mean)
abs_dif_S0_mean <- mean(z_S0_ij_mean)
abs_dif_S0C_mean <- mean(z_S0C_ij_mean)
abs_dif_ground_mean <- mean(z_ground_ij_mean)
abs_dif_mean <- mean(z_ij_mean)
z_sd_S0_mean <- mean(z_S0_ij_sd)
z_sd_S0C_mean <- mean(z_S0C_ij_sd)
z_sd_ground_mean <- mean(z_ground_ij_sd)
z_sd_mean <- mean(z_ij_sd)
length_S0_mean <- mean(l_S0_ij_mean)
length_S0C_mean <- mean(l_S0C_ij_mean)
length_ground_mean <- mean(l_ground_ij_mean)
length_mean <- mean(l_ij_mean)
lengthtrue_S0_mean <- mean(ltrue_S0_ij_mean)
lengthtrue_S0C_mean <- mean(ltrue_S0C_ij_mean)
lengthtrue_ground_mean <- mean(ltrue_ground_ij_mean)
lengthtrue_mean <- mean(ltrue_ij_mean)

#compute sd of asymptotic normality results over S_Omega, S_Omega_C and S_hat_Omega
coverage_rate_S0_sd <- sd(cov_S0_ij_mean)
coverage_rate_S0C_sd <- sd(cov_S0C_ij_mean)
coverage_rate_ground_sd <- sd(cov_ground_ij_mean)
coverage_rate_sd <- sd(cov_ij_mean)
abs_dif_S0_sd <- sd(z_S0_ij_mean)
abs_dif_S0C_sd <- sd(z_S0C_ij_mean)
abs_dif_ground_sd <- sd(z_ground_ij_mean)
abs_dif_sd <- sd(z_ij_mean)
z_sd_S0_sd <- sd(z_S0_ij_sd)
z_sd_S0C_sd <- sd(z_S0C_ij_sd)
z_sd_ground_sd <- sd(z_ground_ij_sd)
z_sd_sd <- sd(z_ij_sd)
length_S0_sd <- sd(l_S0_ij_mean)
length_S0C_sd <- sd(l_S0C_ij_mean)
length_ground_sd <- sd(l_ground_ij_mean)
length_sd <- sd(l_ij_mean)
lengthtrue_S0_sd <- sd(ltrue_S0_ij_mean)
lengthtrue_S0C_sd <- sd(ltrue_S0C_ij_mean)
lengthtrue_ground_sd <- sd(ltrue_ground_ij_mean)
lengthtrue_sd <- sd(ltrue_ij_mean)

rm(z_unbias,l_unbias,ltrue_unbias,cov_unbias)

#Asymptotic normality results of unbias estimator
z_unbias100 <- {}
l_unbias100 <- {}
ltrue_unbias100 <- {}
cov_unbias100 <- {}
for(i in seed){
  dir <- paste("MIO_desparsified_unbias_Gaussian_normality/MIO_desparsified_unbias_Gaussian_band_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
  load(file=dir)
  omega3 <- matrix(1,nrow=p,ncol=p)
  diag(omega3) <- NA
  omega3_sigma_c <- matrix(omega3[which(!is.na(omega3))],nrow=p-1)
  omega3_sigma_c[Ahat] <- NA
  omega3_support<-matrix_generator3(diag(omega3),omega3_sigma_c)
  Ahat1 <- which(is.na(omega3_support))
  #For every replication, find index of support set in intersect of support set of Omega_hat_US under 100 replications
  Ahat2 <- which(Ahat1 %in% ground_support)
  z_unbias100 <- cbind(z_unbias100,z_unbias[Ahat2])
  l_unbias100 <- cbind(l_unbias100,l_unbias[Ahat2])
  ltrue_unbias100 <- cbind(ltrue_unbias100,ltrue_unbias[Ahat2])
  cov_unbias100 <- cbind(cov_unbias100,cov_unbias[Ahat2])
  rm(Ahat,z_unbias,l_unbias,ltrue_unbias,cov_unbias)
}

#For every (i,j), compute mean of asymptotic normality results over 100 replications
z_unbias_ij_mean <- abs(apply(z_unbias100,1,mean))
z_unbias_ij_sd <- apply(z_unbias100,1,sd)
l_unbias_ij_mean <- apply(l_unbias100,1,mean)
ltrue_unbias_ij_mean <- apply(ltrue_unbias100,1,mean)
cov_unbias_ij_mean <- apply(cov_unbias100,1,mean)

#compute mean of asymptotic normality results over S_hat_Omega
coverage_rate_unbias_mean <- mean(cov_unbias_ij_mean)
abs_dif_unbias_mean <- mean(z_unbias_ij_mean)
z_sd_unbias_mean <- mean(z_unbias_ij_sd)
length_unbias_mean <- mean(l_unbias_ij_mean)
lengthtrue_unbias_mean <- mean(ltrue_unbias_ij_mean)

#compute sd of asymptotic normality results over S_hat_Omega
coverage_rate_unbias_sd <- sd(cov_unbias_ij_mean)
abs_dif_unbias_sd <- sd(z_unbias_ij_mean)
z_sd_unbias_sd <- sd(z_unbias_ij_sd)
length_unbias_sd <- sd(l_unbias_ij_mean)
lengthtrue_unbias_sd <- sd(ltrue_unbias_ij_mean)

x1 <- paste(sprintf("%0.3f",length_unbias_mean),"(",sprintf("%0.3f",length_unbias_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_unbias_mean),"(",sprintf("%0.3f",coverage_rate_unbias_sd),")"," & ",
            sprintf("%0.3f",abs_dif_unbias_mean),"(",sprintf("%0.3f",abs_dif_unbias_sd),")"," & ",
            sprintf("%0.3f",z_sd_unbias_mean),"(",sprintf("%0.3f",z_sd_unbias_sd),")"," & ",
            sprintf("%0.3f",lengthtrue_unbias_mean),"(",sprintf("%0.3f",lengthtrue_unbias_sd),")",sep="")
x2 <- paste(sprintf("%0.3f",length_ground_mean),"(",sprintf("%0.3f",length_ground_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_ground_mean),"(",sprintf("%0.3f",coverage_rate_ground_sd),")"," & ",
            sprintf("%0.3f",abs_dif_ground_mean),"(",sprintf("%0.3f",abs_dif_ground_sd),")"," & ",
            sprintf("%0.3f",z_sd_ground_mean),"(",sprintf("%0.3f",z_sd_ground_sd),")"," & ",
            sprintf("%0.3f",lengthtrue_ground_mean),"(",sprintf("%0.3f",lengthtrue_ground_sd),")",sep="")

x3 <- paste(sprintf("%0.3f",length_mean),"(",sprintf("%0.3f",length_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_mean),"(",sprintf("%0.3f",coverage_rate_sd),")"," & ",
            sprintf("%0.3f",abs_dif_mean),"(",sprintf("%0.3f",abs_dif_sd),")"," & ",
            sprintf("%0.3f",z_sd_mean),"(",sprintf("%0.3f",z_sd_sd),")"," & ",
            sprintf("%0.3f",lengthtrue_mean),"(",sprintf("%0.3f",lengthtrue_sd),")",sep="")

x4 <- paste(sprintf("%0.3f",length_S0_mean),"(",sprintf("%0.3f",length_S0_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_S0_mean),"(",sprintf("%0.3f",coverage_rate_S0_sd),")"," & ",
            sprintf("%0.3f",abs_dif_S0_mean),"(",sprintf("%0.3f",abs_dif_S0_sd),")"," & ",
            sprintf("%0.3f",z_sd_S0_mean),"(",sprintf("%0.3f",z_sd_S0_sd),")"," & ",
            sprintf("%0.3f",lengthtrue_S0_mean),"(",sprintf("%0.3f",lengthtrue_S0_sd),")",sep="")

x5 <- paste(sprintf("%0.3f",length_S0C_mean),"(",sprintf("%0.3f",length_S0C_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_S0C_mean),"(",sprintf("%0.3f",coverage_rate_S0C_sd),")"," & ",
            sprintf("%0.3f",abs_dif_S0C_mean),"(",sprintf("%0.3f",abs_dif_S0C_sd),")"," & ",
            sprintf("%0.3f",z_sd_S0C_mean),"(",sprintf("%0.3f",z_sd_S0C_sd),")"," & ",
            sprintf("%0.3f",lengthtrue_S0C_mean),"(",sprintf("%0.3f",lengthtrue_S0C_sd),")",sep="")

x1
x2

x3

x4

x5



