library(Rcpp)
library(RcppArmadillo)
library(huge)
sourceCpp("use_for_all/asdar.cpp")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
seed <- scan("seed.txt",what=integer())
n <- 200 #400
p <- 200 #400,1000,4000

library(BiocParallel)
library(parallel)

ncores=40
mcparam=SnowParam(workers=ncores)

load_function <- function(shu,seed,n,p){
  i <- seed[shu]
  dir <- paste("l1_desparsified_subGaussian_normality/L1_desparsified_subGaussian_band_n",n,"_p",p,"_seed",i,"_normality.RData",sep="")
  load(file=dir)
  
  z_matrix <- matrix(0,nrow=p,ncol=p)
  z_matrix[which(upper.tri(z_matrix, diag = TRUE))] <- z
  z_matrix[lower.tri(z_matrix,diag = FALSE)] <- t(z_matrix)[lower.tri(z_matrix,diag = FALSE)]
  z <- as.vector(z_matrix)
  
  l_matrix <- matrix(0,nrow=p,ncol=p)
  l_matrix[which(upper.tri(l_matrix, diag = TRUE))] <- l
  l_matrix[lower.tri(l_matrix,diag = FALSE)] <- t(l_matrix)[lower.tri(l_matrix,diag = FALSE)]
  l <- as.vector(l_matrix)
  
  cov_matrix <- matrix(0,nrow=p,ncol=p)
  cov_matrix[which(upper.tri(cov_matrix, diag = TRUE))] <- cov
  cov_matrix[lower.tri(cov_matrix,diag = FALSE)] <- t(cov_matrix)[lower.tri(cov_matrix,diag = FALSE)]
  cov <- as.vector(cov_matrix)
  
  return(list(z,l,cov,Ahat))
  
  rm(Ahat,z,l,cov,z_matrix,l_matrix,cov_matrix)
}

iter <- 1:100
result <- bplapply(iter,load_function,seed,n,p)

z_desparsified100 <- matrix(0,nrow=p*p,ncol=100)
l_desparsified100 <- matrix(0,nrow=p*p,ncol=100)
cov_desparsified100 <- matrix(0,nrow=p*p,ncol=100)

for(i in 1:100){
  z_desparsified100[,i] <- result[[i]][[1]]
  l_desparsified100[,i] <- result[[i]][[2]]
  cov_desparsified100[,i] <- result[[i]][[3]]
  print(i)
}

Ahat_ground <- Reduce(intersect, lapply(result[1:100], function(x) x[[4]]))
rm(result)

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
S0 <- which(omega2>0.00001)
#Complementary set of true support set : S_Omega_C
S0C <- which(omega2<=0.00001)
#Intersect of support set of Omega_hat_US under 100 replications : S_hat_Omega
ground_support <- which(is.na(omega_support))

#For every (i,j), compute mean of asymptotic normality results over 100 replications
z_S0_ij_mean <- abs(apply(z_desparsified100[S0,],1,function(x) mean(x, na.rm = TRUE)))
z_S0C_ij_mean <- abs(apply(z_desparsified100[S0C,],1,function(x) mean(x, na.rm = TRUE)))
z_ground_ij_mean <- abs(apply(z_desparsified100[ground_support,],1,function(x) mean(x, na.rm = TRUE)))
z_S0_ij_sd <- apply(z_desparsified100[S0,],1,function(x) sd(x, na.rm = TRUE))
z_S0C_ij_sd <- apply(z_desparsified100[S0C,],1,function(x) sd(x, na.rm = TRUE))
z_ground_ij_sd <- apply(z_desparsified100[ground_support,],1,function(x) sd(x, na.rm = TRUE))
l_S0_ij_mean <- apply(l_desparsified100[S0,],1,function(x) mean(x, na.rm = TRUE))
l_S0C_ij_mean <- apply(l_desparsified100[S0C,],1,function(x) mean(x, na.rm = TRUE))
l_ground_ij_mean <- apply(l_desparsified100[ground_support,],1,function(x) mean(x, na.rm = TRUE))
cov_S0_ij_mean <- apply(cov_desparsified100[S0,],1,function(x) mean(x, na.rm = TRUE))
cov_S0C_ij_mean <- apply(cov_desparsified100[S0C,],1,function(x) mean(x, na.rm = TRUE))
cov_ground_ij_mean <- apply(cov_desparsified100[ground_support,],1,function(x) mean(x, na.rm = TRUE))

#compute mean of asymptotic normality results over S_Omega, S_Omega_C and S_hat_Omega
coverage_rate_S0_mean <- mean(cov_S0_ij_mean,na.rm = TRUE)
coverage_rate_S0C_mean <- mean(cov_S0C_ij_mean,na.rm = TRUE)
coverage_rate_ground_mean <- mean(cov_ground_ij_mean,na.rm = TRUE)
abs_dif_S0_mean <- mean(z_S0_ij_mean,na.rm = TRUE)
abs_dif_S0C_mean <- mean(z_S0C_ij_mean,na.rm = TRUE)
abs_dif_ground_mean <- mean(z_ground_ij_mean,na.rm = TRUE)
z_sd_S0_mean <- mean(z_S0_ij_sd,na.rm = TRUE)
z_sd_S0C_mean <- mean(z_S0C_ij_sd,na.rm = TRUE)
z_sd_ground_mean <- mean(z_ground_ij_sd,na.rm = TRUE)
length_S0_mean <- mean(l_S0_ij_mean,na.rm = TRUE)
length_S0C_mean <- mean(l_S0C_ij_mean,na.rm = TRUE)
length_ground_mean <- mean(l_ground_ij_mean,na.rm = TRUE)

#compute sd of asymptotic normality results over S_Omega, S_Omega_C and S_hat_Omega
coverage_rate_S0_sd <- sd(cov_S0_ij_mean,na.rm = TRUE)
coverage_rate_S0C_sd <- sd(cov_S0C_ij_mean,na.rm = TRUE)
coverage_rate_ground_sd <- sd(cov_ground_ij_mean,na.rm = TRUE)
abs_dif_S0_sd <- sd(z_S0_ij_mean,na.rm = TRUE)
abs_dif_S0C_sd <- sd(z_S0C_ij_mean,na.rm = TRUE)
abs_dif_ground_sd <- sd(z_ground_ij_mean,na.rm = TRUE)
z_sd_S0_sd <- sd(z_S0_ij_sd,na.rm = TRUE)
z_sd_S0C_sd <- sd(z_S0C_ij_sd,na.rm = TRUE)
z_sd_ground_sd <- sd(z_ground_ij_sd,na.rm = TRUE)
length_S0_sd <- sd(l_S0_ij_mean,na.rm = TRUE)
length_S0C_sd <- sd(l_S0C_ij_mean,na.rm = TRUE)
length_ground_sd <- sd(l_ground_ij_mean,na.rm = TRUE)

x1 <- paste(sprintf("%0.3f",length_ground_mean),"(",sprintf("%0.3f",length_ground_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_ground_mean),"(",sprintf("%0.3f",coverage_rate_ground_sd),")"," & ",
            sprintf("%0.3f",abs_dif_ground_mean),"(",sprintf("%0.3f",abs_dif_ground_sd),")"," & ",
            sprintf("%0.3f",z_sd_ground_mean),"(",sprintf("%0.3f",z_sd_ground_sd),")",sep="")

x2 <- paste(sprintf("%0.3f",length_S0_mean),"(",sprintf("%0.3f",length_S0_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_S0_mean),"(",sprintf("%0.3f",coverage_rate_S0_sd),")"," & ",
            sprintf("%0.3f",abs_dif_S0_mean),"(",sprintf("%0.3f",abs_dif_S0_sd),")"," & ",
            sprintf("%0.3f",z_sd_S0_mean),"(",sprintf("%0.3f",z_sd_S0_sd),")",sep="")

x3 <- paste(sprintf("%0.3f",length_S0C_mean),"(",sprintf("%0.3f",length_S0C_sd),")"," & ",
            sprintf("%0.3f",coverage_rate_S0C_mean),"(",sprintf("%0.3f",coverage_rate_S0C_sd),")"," & ",
            sprintf("%0.3f",abs_dif_S0C_mean),"(",sprintf("%0.3f",abs_dif_S0C_sd),")"," & ",
            sprintf("%0.3f",z_sd_S0C_mean),"(",sprintf("%0.3f",z_sd_S0C_sd),")",sep="")

x1

x2

x3

