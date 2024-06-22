#Sys.setlocale('LC_ALL','C')
library(huge)
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
seed <- scan("seed.txt",what=integer())
#Sample size
n <- 200 #400
#Dimension
p <- 200 #400
#Get intersect of support set of Omega_hat_US under 100 replications
Ahat_ground <- c(1:((p-1)*p)) 
#Asymptotic normality results 
z_desparsified100 <- {}
l_desparsified100 <- {}
ltrue_desparsified100 <- {}
cov_desparsified100 <- {}
for(i in seed){
  dir <- paste("l1_desparsified_Gaussian_normality_output/L1_desparsified_Gaussian_band_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
  load(file=dir)
  Ahat_ground <- intersect(Ahat_ground,Ahat)
  z_desparsified100 <- cbind(z_desparsified100,z)
  l_desparsified100 <- cbind(l_desparsified100,l)
  cov_desparsified100 <- cbind(cov_desparsified100,cov)
  rm(Ahat,z,l,cov)
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
#graph= "random"
#set.seed(-1)
#generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, prob=4/p)
#omega2 <- generator$omega
#Sigma2 <- generator$sigma

#Generate Hub graph
#graph= "hub"
#set.seed(-1)
#generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
#omega2 <- generator$omega
#Sigma2 <- generator$sigma

#Generate Cluster graph
#graph= "cluster"
#set.seed(-1)
#generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
#omega2 <- generator$omega
#Sigma2 <- generator$sigma

#True support set : S_Omega
S0 <- which(omega2>0.00001)
#Complementary set of true support set : S_Omega_C
S0C <- which(omega2<=0.00001)
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

x1 <- paste(sprintf("%0.3f",length_ground_mean),"(",sprintf("%0.3f",length_ground_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_ground_mean),"(",sprintf("%0.3f",coverage_rate_ground_sd),")"," & ",
            sprintf("%0.3f",abs_dif_ground_mean),"(",sprintf("%0.3f",abs_dif_ground_sd),")"," & ",
            sprintf("%0.3f",z_sd_ground_mean),"(",sprintf("%0.3f",z_sd_ground_sd),")",sep="")

x2 <- paste(sprintf("%0.3f",length_mean),"(",sprintf("%0.3f",length_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_mean),"(",sprintf("%0.3f",coverage_rate_sd),")"," & ",
            sprintf("%0.3f",abs_dif_mean),"(",sprintf("%0.3f",abs_dif_sd),")"," & ",
            sprintf("%0.3f",z_sd_mean),"(",sprintf("%0.3f",z_sd_sd),")",sep="")

x3 <- paste(sprintf("%0.3f",length_S0_mean),"(",sprintf("%0.3f",length_S0_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_S0_mean),"(",sprintf("%0.3f",coverage_rate_S0_sd),")"," & ",
            sprintf("%0.3f",abs_dif_S0_mean),"(",sprintf("%0.3f",abs_dif_S0_sd),")"," & ",
            sprintf("%0.3f",z_sd_S0_mean),"(",sprintf("%0.3f",z_sd_S0_sd),")",sep="")

x4 <- paste(sprintf("%0.3f",length_S0C_mean),"(",sprintf("%0.3f",length_S0C_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_S0C_mean),"(",sprintf("%0.3f",coverage_rate_S0C_sd),")"," & ",
            sprintf("%0.3f",abs_dif_S0C_mean),"(",sprintf("%0.3f",abs_dif_S0C_sd),")"," & ",
            sprintf("%0.3f",z_sd_S0C_mean),"(",sprintf("%0.3f",z_sd_S0C_sd),")",sep="")

x1

x2

x3

x4

save.image(paste("L1_desparsified_Gaussian_band_n",n,"_p",p,"_normality.RData",sep=""))