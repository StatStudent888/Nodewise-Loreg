library(Rcpp)
library(RcppArmadillo)
library(huge)
sourceCpp("use_for_all/asdar.cpp")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
seed <- scan("seed.txt",what=integer())
n <- 200 #400
p <- 200 #400,1000,4000

load(file=paste(paste("true_length/L0_desparsified_true_Gaussian_subGaussian_band_n",n,"_p",p,"_normality.RData",sep="")))
ltrue_desparsified <- ltrue_Gaussian_desparsified
rm(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified)

# load(file=paste(paste("true_length/L0_desparsified_true_Gaussian_subGaussian_random_n",n,"_p",p,"_normality.RData",sep="")))
# ltrue_desparsified <- ltrue_Gaussian_desparsified
# rm(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified)
# 
# load(file=paste(paste("true_length/L0_desparsified_true_Gaussian_subGaussian_hub_n",n,"_p",p,"_normality.RData",sep="")))
# ltrue_desparsified <- ltrue_Gaussian_desparsified
# rm(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified)
# 
# load(file=paste(paste("true_length/L0_desparsified_true_Gaussian_subGaussian_cluster_n",n,"_p",p,"_normality.RData",sep="")))
# ltrue_desparsified <- ltrue_Gaussian_desparsified
# rm(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified)

library(BiocParallel)
library(parallel)

ncores=40
mcparam=SnowParam(workers=ncores)

load_function_desparsified <- function(shu,seed,n,p){
  i <- seed[shu]
  dir <- paste("l0_desparsified_unbias_Gaussian_normality/L0_desparsified_unbias_Gaussian_band_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
  load(file=dir)
  
  z_matrix <- matrix(0,nrow=p,ncol=p)
  z_matrix[which(upper.tri(z_matrix, diag = TRUE))] <- z_desparsified
  z_matrix[lower.tri(z_matrix,diag = FALSE)] <- t(z_matrix)[lower.tri(z_matrix,diag = FALSE)]
  z_desparsified <- as.vector(z_matrix)
  
  l_matrix <- matrix(0,nrow=p,ncol=p)
  l_matrix[which(upper.tri(l_matrix, diag = TRUE))] <- l_desparsified
  l_matrix[lower.tri(l_matrix,diag = FALSE)] <- t(l_matrix)[lower.tri(l_matrix,diag = FALSE)]
  l_desparsified <- as.vector(l_matrix)
  
  cov_matrix <- matrix(0,nrow=p,ncol=p)
  cov_matrix[which(upper.tri(cov_matrix, diag = TRUE))] <- cov_desparsified
  cov_matrix[lower.tri(cov_matrix,diag = FALSE)] <- t(cov_matrix)[lower.tri(cov_matrix,diag = FALSE)]
  cov_desparsified <- as.vector(cov_matrix)
  
  return(list(z_desparsified,l_desparsified,cov_desparsified,Ahat))
  
  rm(Ahat,z_desparsified,l_desparsified,cov_desparsified,z_matrix,l_matrix,cov_matrix,z_unbias,l_unbias,ltrue_unbias,cov_unbias)
}

# load_function_desparsified <- function(shu,seed,n,p){
#   i <- seed[shu]
#   dir <- paste("l0_desparsified_unbias_Gaussian_normality/L0_desparsified_unbias_Gaussian_random_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
#   load(file=dir)
#   
#   z_matrix <- matrix(0,nrow=p,ncol=p)
#   z_matrix[which(upper.tri(z_matrix, diag = TRUE))] <- z_desparsified
#   z_matrix[lower.tri(z_matrix,diag = FALSE)] <- t(z_matrix)[lower.tri(z_matrix,diag = FALSE)]
#   z_desparsified <- as.vector(z_matrix)
#   
#   l_matrix <- matrix(0,nrow=p,ncol=p)
#   l_matrix[which(upper.tri(l_matrix, diag = TRUE))] <- l_desparsified
#   l_matrix[lower.tri(l_matrix,diag = FALSE)] <- t(l_matrix)[lower.tri(l_matrix,diag = FALSE)]
#   l_desparsified <- as.vector(l_matrix)
#   
#   cov_matrix <- matrix(0,nrow=p,ncol=p)
#   cov_matrix[which(upper.tri(cov_matrix, diag = TRUE))] <- cov_desparsified
#   cov_matrix[lower.tri(cov_matrix,diag = FALSE)] <- t(cov_matrix)[lower.tri(cov_matrix,diag = FALSE)]
#   cov_desparsified <- as.vector(cov_matrix)
#   
#   return(list(z_desparsified,l_desparsified,cov_desparsified,Ahat))
#   
#   rm(Ahat,z_desparsified,l_desparsified,cov_desparsified,z_matrix,l_matrix,cov_matrix,z_unbias,l_unbias,ltrue_unbias,cov_unbias)
# }
# 
# load_function_desparsified <- function(shu,seed,n,p){
#   i <- seed[shu]
#   dir <- paste("l0_desparsified_unbias_Gaussian_normality/L0_desparsified_unbias_Gaussian_hub_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
#   load(file=dir)
#   
#   z_matrix <- matrix(0,nrow=p,ncol=p)
#   z_matrix[which(upper.tri(z_matrix, diag = TRUE))] <- z_desparsified
#   z_matrix[lower.tri(z_matrix,diag = FALSE)] <- t(z_matrix)[lower.tri(z_matrix,diag = FALSE)]
#   z_desparsified <- as.vector(z_matrix)
#   
#   l_matrix <- matrix(0,nrow=p,ncol=p)
#   l_matrix[which(upper.tri(l_matrix, diag = TRUE))] <- l_desparsified
#   l_matrix[lower.tri(l_matrix,diag = FALSE)] <- t(l_matrix)[lower.tri(l_matrix,diag = FALSE)]
#   l_desparsified <- as.vector(l_matrix)
#   
#   cov_matrix <- matrix(0,nrow=p,ncol=p)
#   cov_matrix[which(upper.tri(cov_matrix, diag = TRUE))] <- cov_desparsified
#   cov_matrix[lower.tri(cov_matrix,diag = FALSE)] <- t(cov_matrix)[lower.tri(cov_matrix,diag = FALSE)]
#   cov_desparsified <- as.vector(cov_matrix)
#   
#   return(list(z_desparsified,l_desparsified,cov_desparsified,Ahat))
#   
#   rm(Ahat,z_desparsified,l_desparsified,cov_desparsified,z_matrix,l_matrix,cov_matrix,z_unbias,l_unbias,ltrue_unbias,cov_unbias)
# }
# 
# load_function_desparsified <- function(shu,seed,n,p){
#   i <- seed[shu]
#   dir <- paste("l0_desparsified_unbias_Gaussian_normality/L0_desparsified_unbias_Gaussian_cluster_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
#   load(file=dir)
#   
#   z_matrix <- matrix(0,nrow=p,ncol=p)
#   z_matrix[which(upper.tri(z_matrix, diag = TRUE))] <- z_desparsified
#   z_matrix[lower.tri(z_matrix,diag = FALSE)] <- t(z_matrix)[lower.tri(z_matrix,diag = FALSE)]
#   z_desparsified <- as.vector(z_matrix)
#   
#   l_matrix <- matrix(0,nrow=p,ncol=p)
#   l_matrix[which(upper.tri(l_matrix, diag = TRUE))] <- l_desparsified
#   l_matrix[lower.tri(l_matrix,diag = FALSE)] <- t(l_matrix)[lower.tri(l_matrix,diag = FALSE)]
#   l_desparsified <- as.vector(l_matrix)
#   
#   cov_matrix <- matrix(0,nrow=p,ncol=p)
#   cov_matrix[which(upper.tri(cov_matrix, diag = TRUE))] <- cov_desparsified
#   cov_matrix[lower.tri(cov_matrix,diag = FALSE)] <- t(cov_matrix)[lower.tri(cov_matrix,diag = FALSE)]
#   cov_desparsified <- as.vector(cov_matrix)
#   
#   return(list(z_desparsified,l_desparsified,cov_desparsified,Ahat))
#   
#   rm(Ahat,z_desparsified,l_desparsified,cov_desparsified,z_matrix,l_matrix,cov_matrix,z_unbias,l_unbias,ltrue_unbias,cov_unbias)
# }

iter <- 1:100
result_desparsified <- bplapply(iter,load_function_desparsified,seed,n,p)

z_desparsified100 <- matrix(0,nrow=p*p,ncol=100)
l_desparsified100 <- matrix(0,nrow=p*p,ncol=100)
ltrue_desparsified100 <- matrix(0,nrow=p*p,ncol=100)
cov_desparsified100 <- matrix(0,nrow=p*p,ncol=100)

for(i in 1:100){
  z_desparsified100[,i] <- result_desparsified[[i]][[1]]
  l_desparsified100[,i] <- result_desparsified[[i]][[2]]
  cov_desparsified100[,i] <- result_desparsified[[i]][[3]]
  ltrue_desparsified100[,i] <- ltrue_desparsified
  # print(i)
}

Ahat_ground <- Reduce(intersect, lapply(result_desparsified[1:100], function(x) x[[4]]))
rm(result_desparsified)

omega <- matrix(1,nrow=p,ncol=p)
diag(omega) <- NA
omega_sigma_c <- matrix(omega[which(!is.na(omega))],nrow=p-1)
omega_sigma_c[Ahat_ground] <- NA
omega_support<-matrix_generator3(diag(omega),omega_sigma_c)

#Generate Random graph
set.seed(-1)
omega2 <- omega_generator1(p,1,0.5,0.3)
if(p %in% c(200,400)){
  Sigma2 <- solve(omega2)
} else {Sigma2 <- matrix_inverse(omega2)}

#Generate Band graph
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
# Sigma2<-generator$sigma
# omega2<-generator$omega

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
z_S0_ij_sd <- apply(z_desparsified100[S0,],1,sd)
z_S0C_ij_sd <- apply(z_desparsified100[S0C,],1,sd)
z_ground_ij_sd <- apply(z_desparsified100[ground_support,],1,sd)
l_S0_ij_mean <- apply(l_desparsified100[S0,],1,mean)
l_S0C_ij_mean <- apply(l_desparsified100[S0C,],1,mean)
l_ground_ij_mean <- apply(l_desparsified100[ground_support,],1,mean)
ltrue_S0_ij_mean <- apply(ltrue_desparsified100[S0,],1,mean)
ltrue_S0C_ij_mean <- apply(ltrue_desparsified100[S0C,],1,mean)
ltrue_ground_ij_mean <- apply(ltrue_desparsified100[ground_support,],1,mean)
cov_S0_ij_mean <- apply(cov_desparsified100[S0,],1,mean)
cov_S0C_ij_mean <- apply(cov_desparsified100[S0C,],1,mean)
cov_ground_ij_mean <- apply(cov_desparsified100[ground_support,],1,mean)

#compute mean of asymptotic normality results over S_Omega, S_Omega_C and S_hat_Omega
coverage_rate_S0_mean <- mean(cov_S0_ij_mean)
coverage_rate_S0C_mean <- mean(cov_S0C_ij_mean)
coverage_rate_ground_mean <- mean(cov_ground_ij_mean)
abs_dif_S0_mean <- mean(z_S0_ij_mean)
abs_dif_S0C_mean <- mean(z_S0C_ij_mean)
abs_dif_ground_mean <- mean(z_ground_ij_mean)
z_sd_S0_mean <- mean(z_S0_ij_sd)
z_sd_S0C_mean <- mean(z_S0C_ij_sd)
z_sd_ground_mean <- mean(z_ground_ij_sd)
length_S0_mean <- mean(l_S0_ij_mean)
length_S0C_mean <- mean(l_S0C_ij_mean)
length_ground_mean <- mean(l_ground_ij_mean)
lengthtrue_S0_mean <- mean(ltrue_S0_ij_mean)
lengthtrue_S0C_mean <- mean(ltrue_S0C_ij_mean)
lengthtrue_ground_mean <- mean(ltrue_ground_ij_mean)

#compute sd of asymptotic normality results over S_Omega, S_Omega_C and S_hat_Omega
coverage_rate_S0_sd <- sd(cov_S0_ij_mean)
coverage_rate_S0C_sd <- sd(cov_S0C_ij_mean)
coverage_rate_ground_sd <- sd(cov_ground_ij_mean)
abs_dif_S0_sd <- sd(z_S0_ij_mean)
abs_dif_S0C_sd <- sd(z_S0C_ij_mean)
abs_dif_ground_sd <- sd(z_ground_ij_mean)
z_sd_S0_sd <- sd(z_S0_ij_sd)
z_sd_S0C_sd <- sd(z_S0C_ij_sd)
z_sd_ground_sd <- sd(z_ground_ij_sd)
length_S0_sd <- sd(l_S0_ij_mean)
length_S0C_sd <- sd(l_S0C_ij_mean)
length_ground_sd <- sd(l_ground_ij_mean)
lengthtrue_S0_sd <- sd(ltrue_S0_ij_mean)
lengthtrue_S0C_sd <- sd(ltrue_S0C_ij_mean)
lengthtrue_ground_sd <- sd(ltrue_ground_ij_mean)

rm(z_desparsified100,l_desparsified100,cov_desparsified100,Ahat_ground)

load_function_unbias <- function(shu,seed,n,p){
  i <- seed[shu]
  dir <- paste("l0_desparsified_unbias_Gaussian_normality/L0_desparsified_unbias_Gaussian_band_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
  load(file=dir)
  omega3 <- matrix(1,nrow=p,ncol=p)
  diag(omega3) <- NA
  omega3_sigma_c <- matrix(omega3[which(!is.na(omega3))],nrow=p-1)
  omega3_sigma_c[Ahat] <- NA
  omega3_support<-matrix_generator3(diag(omega3),omega3_sigma_c)
  Ahat1 <- which(is.na(omega3_support))
  #For every replication, find index of support set in intersect of support set of Omega_hat_US under 100 replications
  Ahat2 <- which(Ahat1 %in% ground_support)
  
  return(list(z_unbias=z_unbias[Ahat2],l_unbias=l_unbias[Ahat2],ltrue_unbias=ltrue_unbias[Ahat2],cov_unbias=cov_unbias[Ahat2]))
  
  rm(Ahat,z_desparsified,l_desparsified,cov_desparsified,z_unbias,l_unbias,ltrue_unbias,cov_unbias)
}

# load_function_unbias <- function(shu,seed,n,p){
#   i <- seed[shu]
#   dir <- paste("l0_desparsified_unbias_Gaussian_normality/L0_desparsified_unbias_Gaussian_random_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
#   load(file=dir)
#   omega3 <- matrix(1,nrow=p,ncol=p)
#   diag(omega3) <- NA
#   omega3_sigma_c <- matrix(omega3[which(!is.na(omega3))],nrow=p-1)
#   omega3_sigma_c[Ahat] <- NA
#   omega3_support<-matrix_generator3(diag(omega3),omega3_sigma_c)
#   Ahat1 <- which(is.na(omega3_support))
#   #For every replication, find index of support set in intersect of support set of Omega_hat_US under 100 replications
#   Ahat2 <- which(Ahat1 %in% ground_support)
#   
#   return(list(z_unbias=z_unbias[Ahat2],l_unbias=l_unbias[Ahat2],ltrue_unbias=ltrue_unbias[Ahat2],cov_unbias=cov_unbias[Ahat2]))
#   
#   rm(Ahat,z_desparsified,l_desparsified,cov_desparsified,z_unbias,l_unbias,ltrue_unbias,cov_unbias)
# }
# 
# load_function_unbias <- function(shu,seed,n,p){
#   i <- seed[shu]
#   dir <- paste("l0_desparsified_unbias_Gaussian_normality/L0_desparsified_unbias_Gaussian_hub_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
#   load(file=dir)
#   omega3 <- matrix(1,nrow=p,ncol=p)
#   diag(omega3) <- NA
#   omega3_sigma_c <- matrix(omega3[which(!is.na(omega3))],nrow=p-1)
#   omega3_sigma_c[Ahat] <- NA
#   omega3_support<-matrix_generator3(diag(omega3),omega3_sigma_c)
#   Ahat1 <- which(is.na(omega3_support))
#   #For every replication, find index of support set in intersect of support set of Omega_hat_US under 100 replications
#   Ahat2 <- which(Ahat1 %in% ground_support)
#   
#   return(list(z_unbias=z_unbias[Ahat2],l_unbias=l_unbias[Ahat2],ltrue_unbias=ltrue_unbias[Ahat2],cov_unbias=cov_unbias[Ahat2]))
#   
#   rm(Ahat,z_desparsified,l_desparsified,cov_desparsified,z_unbias,l_unbias,ltrue_unbias,cov_unbias)
# }
# 
# load_function_unbias <- function(shu,seed,n,p){
#   i <- seed[shu]
#   dir <- paste("l0_desparsified_unbias_Gaussian_normality/L0_desparsified_unbias_Gaussian_cluster_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
#   load(file=dir)
#   omega3 <- matrix(1,nrow=p,ncol=p)
#   diag(omega3) <- NA
#   omega3_sigma_c <- matrix(omega3[which(!is.na(omega3))],nrow=p-1)
#   omega3_sigma_c[Ahat] <- NA
#   omega3_support<-matrix_generator3(diag(omega3),omega3_sigma_c)
#   Ahat1 <- which(is.na(omega3_support))
#   #For every replication, find index of support set in intersect of support set of Omega_hat_US under 100 replications
#   Ahat2 <- which(Ahat1 %in% ground_support)
#   
#   return(list(z_unbias=z_unbias[Ahat2],l_unbias=l_unbias[Ahat2],ltrue_unbias=ltrue_unbias[Ahat2],cov_unbias=cov_unbias[Ahat2]))
#   
#   rm(Ahat,z_desparsified,l_desparsified,cov_desparsified,z_unbias,l_unbias,ltrue_unbias,cov_unbias)
# }

result_unbias <- bplapply(iter,load_function_unbias,seed,n,p)

z_unbias100 <- {}
l_unbias100 <- {}
ltrue_unbias100 <- {}
cov_unbias100 <- {}

for(i in 1:100){
  z_unbias <- result_unbias[[i]][[1]]
  l_unbias <- result_unbias[[i]][[2]]
  ltrue_unbias <- result_unbias[[i]][[3]]
  cov_unbias <- result_unbias[[i]][[4]]
  
  z_unbias100 <- cbind(z_unbias100,z_unbias)
  l_unbias100 <- cbind(l_unbias100,l_unbias)
  ltrue_unbias100 <- cbind(ltrue_unbias100,ltrue_unbias)
  cov_unbias100 <- cbind(cov_unbias100,cov_unbias)
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

x1 <- paste(sprintf("%0.3f",lengthtrue_unbias_mean),"(",sprintf("%0.3f",lengthtrue_unbias_sd),")"," & ",
            sprintf("%0.3f",length_unbias_mean),"(",sprintf("%0.3f",length_unbias_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_unbias_mean),"(",sprintf("%0.3f",coverage_rate_unbias_sd),")"," & ",
            sprintf("%0.3f",abs_dif_unbias_mean),"(",sprintf("%0.3f",abs_dif_unbias_sd),")"," & ",
            sprintf("%0.3f",z_sd_unbias_mean),"(",sprintf("%0.3f",z_sd_unbias_sd),")",sep="")

x2 <- paste(sprintf("%0.3f",lengthtrue_ground_mean),"(",sprintf("%0.3f",lengthtrue_ground_sd),")"," & ",
            sprintf("%0.3f",length_ground_mean),"(",sprintf("%0.3f",length_ground_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_ground_mean),"(",sprintf("%0.3f",coverage_rate_ground_sd),")"," & ",
            sprintf("%0.3f",abs_dif_ground_mean),"(",sprintf("%0.3f",abs_dif_ground_sd),")"," & ",
            sprintf("%0.3f",z_sd_ground_mean),"(",sprintf("%0.3f",z_sd_ground_sd),")",sep="")

x3 <- paste(sprintf("%0.3f",lengthtrue_S0_mean),"(",sprintf("%0.3f",lengthtrue_S0_sd),")"," & ",
            sprintf("%0.3f",length_S0_mean),"(",sprintf("%0.3f",length_S0_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_S0_mean),"(",sprintf("%0.3f",coverage_rate_S0_sd),")"," & ",
            sprintf("%0.3f",abs_dif_S0_mean),"(",sprintf("%0.3f",abs_dif_S0_sd),")"," & ",
            sprintf("%0.3f",z_sd_S0_mean),"(",sprintf("%0.3f",z_sd_S0_sd),")",sep="")

x4 <- paste(sprintf("%0.3f",lengthtrue_S0C_mean),"(",sprintf("%0.3f",lengthtrue_S0C_sd),")"," & ",
            sprintf("%0.3f",length_S0C_mean),"(",sprintf("%0.3f",length_S0C_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_S0C_mean),"(",sprintf("%0.3f",coverage_rate_S0C_sd),")"," & ",
            sprintf("%0.3f",abs_dif_S0C_mean),"(",sprintf("%0.3f",abs_dif_S0C_sd),")"," & ",
            sprintf("%0.3f",z_sd_S0C_mean),"(",sprintf("%0.3f",z_sd_S0C_sd),")",sep="")


x1
x2

x3

x4


