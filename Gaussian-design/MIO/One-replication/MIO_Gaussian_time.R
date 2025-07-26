library(reticulate)
conda_list()
use_condaenv("base", required = TRUE)
py_config()
library(Rcpp)
library(RcppArmadillo)
library(huge)
library(MASS)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
py_run_string("import sys; sys.path.append('use_for_all/subset_selection_with_shrinkage-master/python/algorithms')")
source_python("use_for_all/subset_selection_with_shrinkage-master/python/algorithms/MIO.py")
source("use_for_all/amio-end.R")
sourceCpp("use_for_all/asdar.cpp")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
source("use_for_all/evaluation.R")
source("use_for_all/useful_function.R")

#Sample size
n=200 
#Dimension
p=200 

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

SEED<-scan("seed.txt")
seed = SEED[1] #1-100
set.seed(seed)
#Generate Gaussian data
data <- mvrnorm(n,Mean,Sigma2)

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
    best_ebeta<-amio(X,y,20,n,p,'l1',5)
    sigmak<-n/(base::norm(y-X%*%best_ebeta,"2")^2)
    sigma_ck<-(-1)*sigmak*best_ebeta  
    sigma<- cbind(sigma,sigmak)
    sigma_c<- cbind(sigma_c,sigma_ck)
  }
  sigma<-sigma[-1]
  sigma_c<-sigma_c[,-1]
  MIO_F_sigma1<-matrix_generator1(sigma,sigma_c) }) 

#Support recovery performance
eva <- evaluation(omega2,MIO_F_sigma1)

precision <- eva$precision
sensitivity <- eva$sensitivity
specificity <- eva$specificity
MCC <- eva$MCC

#Matrix norm loss
MIO_norm1<-base::norm(MIO_F_sigma1-omega2,"1")
MIO_norm2<-base::norm(MIO_F_sigma1-omega2,"2")
MIO_normF<-base::norm(MIO_F_sigma1-omega2,"F")
MIO_normmax<-base::norm(MIO_F_sigma1-omega2,"M")

#Running time
MIO_time<-timek

save(precision,sensitivity,specificity,MCC,MIO_norm1,MIO_norm2,MIO_normF,MIO_normmax,MIO_time,file=paste("MIO_Gaussian_band_n",n,"_p",p,"_seed",seed,".RData",sep=""))

save(precision,sensitivity,specificity,MCC,MIO_norm1,MIO_norm2,MIO_normF,MIO_normmax,MIO_time,file=paste("MIO_Gaussian_random_n",n,"_p",p,"_seed",seed,".RData",sep=""))

save(precision,sensitivity,specificity,MCC,MIO_norm1,MIO_norm2,MIO_normF,MIO_normmax,MIO_time,file=paste("MIO_Gaussian_hub_n",n,"_p",p,"_seed",seed,".RData",sep=""))

save(precision,sensitivity,specificity,MCC,MIO_norm1,MIO_norm2,MIO_normF,MIO_normmax,MIO_time,file=paste("MIO_Gaussian_cluster_n",n,"_p",p,"_seed",seed,".RData",sep=""))
