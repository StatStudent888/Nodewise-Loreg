library(huge)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
sourceCpp("use_for_all/compute_desparsified.cpp")

COL <- 1:4000
column <- COL[1] #1-4000

#Sample size
n=200 #400
#Dimension
p=4000 
#FDR level
alpha <- 0.05

# #Generate Band graph
set.seed(-1)
omega2 <- omega_generator1(p,1,0.5,0.3)
Sigma2 <- matrix_inverse(omega2)
Sigma_sqrt <- matrix_sqrt(Sigma2)

#Generate Random graph
# graph= "random"
# set.seed(-1)
# generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, prob=4/p)
# omega2 <- generator$omega
# Sigma2 <- generator$sigma
# Sigma_sqrt <- matrix_sqrt(Sigma2)

#Generate Hub graph
# graph= "hub"
# set.seed(-1)
# generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
# omega2 <- generator$omega
# Sigma2 <- generator$sigma
# Sigma_sqrt <- matrix_sqrt(Sigma2)

#Generate Cluster graph
# graph= "cluster"
# set.seed(-1)
# generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
# omega2 <- generator$omega
# Sigma2 <- generator$sigma
# Sigma_sqrt <- matrix_sqrt(Sigma2)

result <- compute_desparsified(column-1,omega2,Sigma_sqrt,alpha,n,p)

ltrue_Gaussian_desparsified <- result$ltrue_Gaussian_desparsified
ltrue_subGaussian_desparsified <- result$ltrue_subGaussian_desparsified

save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_Gaussian_subGaussian_band_n",n,"_p",p,"_col",column,"_true_length.RData",sep=""))

# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_Gaussian_subGaussian_random_n",n,"_p",p,"_col",column,"_true_length.RData",sep=""))
# 
# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_Gaussian_subGaussian_hub_n",n,"_p",p,"_col",column,"_true_length.RData",sep=""))
# 
# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_Gaussian_subGaussian_cluster_n",n,"_p",p,"_col",column,"_true_length.RData",sep=""))
