library(reticulate)
conda_list()
use_condaenv("base", required = TRUE)
py_config()
library(Rcpp)
library(RcppArmadillo)
library(huge)
library(MASS)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
source_python("use_for_all/l0bnb-master/l0bnb/regpath.py")
source("use_for_all/abnb-end.R")
sourceCpp("use_for_all/asdar.cpp")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
source("use_for_all/evaluation.R")
source("use_for_all/useful_function.R")

#Sample size
n=200 
#Dimension
p=200 
#FDR level
alpha <- 0.05

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

Mean <- rep(0,p)

#Generate Sub-Gaussian data
Sigma_sqrt <- matrix_sqrt(Sigma2)
U <- matrix(0,n,p)

SEED<-scan("seed.txt")
seed = SEED[1] #1-100
set.seed(seed)

for(h in 1:n){
  for(q in 1:p){
    U[h,q] <- runif(1,-sqrt(3),sqrt(3))
  }
}
data <- matrixMultiply(U,Sigma_sqrt)

#Estimate precision matrix and compute running time
j<-0
sigma<-matrix()
sigma_c<-matrix(1,p-1,1)
timek <- system.time({
  while(j<p){
    j<-j+1
    X<-data[,-j]
    y<-data[,j]
    X<-as.matrix(X) 
    y<-as.vector(y)
    best_ebeta<-abnb(X,y,20,n,p,m_multi=1.2,time_lim=5,tol=0.01,lambda0r=0.9,max_iter=10)
    sigmak<-n/(base::norm(y-X%*%best_ebeta,"2")^2)
    sigma_ck<-(-1)*sigmak*best_ebeta  
    sigma<- cbind(sigma,sigmak)
    sigma_c<- cbind(sigma_c,sigma_ck)
  }
  sigma<-sigma[-1]
  sigma_c<-sigma_c[,-1]
  BnB_F_sigma1<-matrix_generator1(sigma,sigma_c) }) 

#Support recovery performance
eva <- evaluation(omega2,BnB_F_sigma1)

precision <- eva$precision
sensitivity <- eva$sensitivity
specificity <- eva$specificity
MCC <- eva$MCC

#Matrix norm loss
BnB_norm1<-base::norm(BnB_F_sigma1-omega2,"1")
BnB_norm2<-base::norm(BnB_F_sigma1-omega2,"2")
BnB_normF<-base::norm(BnB_F_sigma1-omega2,"F")
BnB_normmax<-base::norm(BnB_F_sigma1-omega2,"M")

#Running time
BnB_time<-timek

save(precision,sensitivity,specificity,MCC,BnB_norm1,BnB_norm2,BnB_normF,BnB_normmax,BnB_time,file=paste("BnB_subGaussian_band_n",n,"_p",p,"_seed",seed,".RData",sep=""))

# save(precision,sensitivity,specificity,MCC,BnB_norm1,BnB_norm2,BnB_normF,BnB_normmax,BnB_time,file=paste("BnB_subGaussian_random_n",n,"_p",p,"_seed",seed,".RData",sep=""))
# 
# save(precision,sensitivity,specificity,MCC,BnB_norm1,BnB_norm2,BnB_normF,BnB_normmax,BnB_time,file=paste("BnB_subGaussian_hub_n",n,"_p",p,"_seed",seed,".RData",sep=""))
# 
# save(precision,sensitivity,specificity,MCC,BnB_norm1,BnB_norm2,BnB_normF,BnB_normmax,BnB_time,file=paste("BnB_subGaussian_cluster_n",n,"_p",p,"_seed",seed,".RData",sep=""))






