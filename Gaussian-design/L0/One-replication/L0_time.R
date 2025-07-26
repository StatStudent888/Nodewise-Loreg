library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(huge)
library(MASS)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("use_for_all/asdar_youhua_Eigen.cpp")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
source("use_for_all/useful_function.R")

#Sample size
n=200 #400
#Dimension
p=200 #400,1000,4000

#Generate Band graph
set.seed(-1)
omega2 <- omega_generator1(p,1,0.5,0.3)
if(p %in% c(200,400)){
  Sigma2 <- solve(omega2)
} else {Sigma2 <- matrix_inverse(omega2)}

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

Mean <- rep(0,p)

SEED<-scan("seed.txt")
seed = SEED[1] #1-100
set.seed(seed)

#Generate Gaussian data
data_Gaussian <- mvrnorm(n,Mean,Sigma2)
result_Gaussian <- normalizeMatrix(data_Gaussian)
data_Gaussian_normalize <- result_Gaussian$normalizedX
diag_Gaussian <- result_Gaussian$norms

#Generate Sub-Gaussian data
Sigma_sqrt <- matrix_sqrt(Sigma2)
U <- matrix(0,n,p)
set.seed(seed)
for(h in 1:n){
  for(q in 1:p){
    U[h,q] <- runif(1,-sqrt(3),sqrt(3))
  }
}
data_subGaussian <- U%*%Sigma_sqrt
result_subGaussian <- normalizeMatrix(data_subGaussian)
data_subGaussian_normalize <- result_subGaussian$normalizedX
diag_subGaussian <- result_subGaussian$norms

#Compute running time
timek_Gaussian<-system.time({
  result <- process(data_Gaussian,data_Gaussian_normalize,diag_Gaussian,n,p)
  sdar_F_sigma1_Gaussian<-matrix_generator1(as.vector(result$sigma),result$sigma_c) }) 

timek_subGaussian<-system.time({
  result <- process(data_subGaussian,data_subGaussian_normalize,diag_subGaussian,n,p)
  sdar_F_sigma1_subGaussian<-matrix_generator1(as.vector(result$sigma),result$sigma_c) }) 


save(timek_Gaussian,timek_subGaussian,file=paste("L0_band_n",n,"_p",p,"_seed",seed,"_time.RData",sep=""))

# save(timek_Gaussian,timek_subGaussian,file=paste("L0_random_n",n,"_p",p,"_seed",seed,"_time.RData",sep=""))
# 
# save(timek_Gaussian,timek_subGaussian,file=paste("L0_hub_n",n,"_p",p,"_seed",seed,"_time.RData",sep=""))
# 
# save(timek_Gaussian,timek_subGaussian,file=paste("L0_cluster_n",n,"_p",p,"_seed",seed,"_time.RData",sep=""))



