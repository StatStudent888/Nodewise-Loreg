library(Rcpp)
library(RcppArmadillo)
library(huge)
library(MASS)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("use_for_all/asdar.cpp")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
source("use_for_all/evaluation.R")
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
data <- mvrnorm(n,Mean,Sigma2)
X_plus <- as.matrix(data)
X_plus2 <- matrixMultiply(t(X_plus),X_plus)/n
  
result <- normalizeMatrix(data)
data_normalize <- result$normalizedX
diag <- result$norms

#Estimate precision matrix and compute the number of iteration
j<-0
iter_count <- {}
iter_flag <- {}
while(j<p){
  j<-j+1
  X<-data[,-j]
  y<-data[,j]
  nX <- data_normalize[,-j]
  dx <- diag[-j]
  y<-as.vector(y)
  result_list <- asdar_iter(X,nX,dx,y,20,n,p)
  countk <- result_list$iter_count
  flagk <- result_list$iter_flag
  iter_count <- c(iter_count,countk)
  iter_flag <- c(iter_flag,flagk)
}
count_mean <- mean(iter_count)
count_sd <- sd(iter_count)
count_median <- median(iter_count)
flag_mean <- mean(iter_flag) 

save(count_mean,count_sd,count_median,flag_mean,iter_count,file=paste("L0_Gaussian_band_n",n,"_p",p,"_seed",seed,"_NumOfIter.RData",sep=""))

# save(count_mean,count_sd,count_median,flag_mean,iter_count,file=paste("L0_Gaussian_random_n",n,"_p",p,"_seed",seed,"_NumOfIter.RData",sep=""))
# 
# save(count_mean,count_sd,count_median,flag_mean,iter_count,file=paste("L0_Gaussian_hub_n",n,"_p",p,"_seed",seed,"_NumOfIter.RData",sep=""))
# 
# save(count_mean,count_sd,count_median,flag_mean,iter_count,file=paste("L0_Gaussian_cluster_n",n,"_p",p,"_seed",seed,"_NumOfIter.RData",sep=""))

