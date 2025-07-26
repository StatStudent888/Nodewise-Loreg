library(huge)
library(Rcpp)
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
library(BiocParallel)
library(parallel)
sourceCpp("use_for_all/compute_desparsified.cpp")

#Sample size
n=200 #400
#Dimension
p=200 #400,1000
#FDR level
alpha <- 0.05

ncores = 40
mcparam = SnowParam(workers = ncores)

# #Generate Band graph
set.seed(-1)
omega2 <- omega_generator1(p,1,0.5,0.3)
if(p %in% c(200,400)){
  Sigma2 <- solve(omega2)
} else {Sigma2 <- matrix_inverse(omega2)}
Sigma_sqrt <- matrix_sqrt(Sigma2)

#Generate Random graph
#graph= "random"
#set.seed(-1)
#generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, prob=4/p)
#omega2 <- generator$omega
#Sigma2 <- generator$sigma
#Sigma_sqrt <- matrix_sqrt(Sigma2)

#Generate Hub graph
#graph= "hub"
#set.seed(-1)
#generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
#omega2 <- generator$omega
#Sigma2 <- generator$sigma
#Sigma_sqrt <- matrix_sqrt(Sigma2)

#Generate Cluster graph
# graph= "cluster"
# set.seed(-1)
# generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
# omega2 <- generator$omega
# Sigma2 <- generator$sigma
# Sigma_sqrt <- matrix_sqrt(Sigma2)

S0 <- which(abs(omega2)>0.00001)
S0C <- which(abs(omega2)<=0.00001)

desparsified <- function(j,omega2,Sigma_sqrt,alpha,n,p){
  library(Rcpp)
  sourceCpp("use_for_all/compute_desparsified.cpp")
  result <- compute_desparsified(j,omega2,Sigma_sqrt,alpha,n,p)
  return(list(result$ltrue_Gaussian_desparsified,result$ltrue_subGaussian_desparsified))
}

results <- bplapply(c(0:(p-1)),desparsified,omega2,Sigma_sqrt,alpha,n,p,BPPARAM = mcparam)

ltrue_Gaussian_desparsified <- {}
ltrue_subGaussian_desparsified <- {}
for(i in 1:p){
  ltrue_Gaussian_desparsified <- c(ltrue_Gaussian_desparsified,results[[i]][[1]])
  ltrue_subGaussian_desparsified <- c(ltrue_subGaussian_desparsified,results[[i]][[2]])
}

Gaussian_matrix <- matrix(0,nrow=p,ncol=p)
Gaussian_matrix[which(upper.tri(Gaussian_matrix, diag = TRUE))] <- ltrue_Gaussian_desparsified
Gaussian_matrix[lower.tri(Gaussian_matrix,diag = FALSE)] <- t(Gaussian_matrix)[lower.tri(Gaussian_matrix,diag = FALSE)]
ltrue_Gaussian_desparsified <- as.vector(Gaussian_matrix)

subGaussian_matrix <- matrix(0,nrow=p,ncol=p)
subGaussian_matrix[which(upper.tri(subGaussian_matrix, diag = TRUE))] <- ltrue_subGaussian_desparsified
subGaussian_matrix[lower.tri(subGaussian_matrix,diag = FALSE)] <- t(subGaussian_matrix)[lower.tri(subGaussian_matrix,diag = FALSE)]
ltrue_subGaussian_desparsified <- as.vector(subGaussian_matrix)


lengthtrue_Gaussian_S0_mean <- mean(ltrue_Gaussian_desparsified[S0])
lengthtrue_Gaussian_S0C_mean <- mean(ltrue_Gaussian_desparsified[S0C])

lengthtrue_Gaussian_S0_sd <- sd(ltrue_Gaussian_desparsified[S0])
lengthtrue_Gaussian_S0C_sd <- sd(ltrue_Gaussian_desparsified[S0C])

lengthtrue_subGaussian_S0_mean <- mean(ltrue_subGaussian_desparsified[S0])
lengthtrue_subGaussian_S0C_mean <- mean(ltrue_subGaussian_desparsified[S0C])

lengthtrue_subGaussian_S0_sd <- sd(ltrue_subGaussian_desparsified[S0])
lengthtrue_subGaussian_S0C_sd <- sd(ltrue_subGaussian_desparsified[S0C])

x1 <- paste(sprintf("%0.3f",lengthtrue_Gaussian_S0_mean),"(",sprintf("%0.3f",lengthtrue_Gaussian_S0_sd),")",sep="")

x2 <- paste(sprintf("%0.3f",lengthtrue_Gaussian_S0C_mean),"(",sprintf("%0.3f",lengthtrue_Gaussian_S0C_sd),")",sep="")

x3 <- paste(sprintf("%0.3f",lengthtrue_subGaussian_S0_mean),"(",sprintf("%0.3f",lengthtrue_subGaussian_S0_sd),")",sep="")

x4 <- paste(sprintf("%0.3f",lengthtrue_subGaussian_S0C_mean),"(",sprintf("%0.3f",lengthtrue_subGaussian_S0C_sd),")",sep="")

x1

x2

x3

x4

save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_band_n",n,"_p",p,"_normality.RData",sep=""))

# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_random_n",n,"_p",p,"_normality.RData",sep=""))
# 
# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_hub_n",n,"_p",p,"_normality.RData",sep=""))
# 
# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_cluster_n",n,"_p",p,"_normality.RData",sep=""))


