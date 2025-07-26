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
source("use_for_all/FDR-R-Code/epsest.func.R")
source("use_for_all/FDR-R-Code/lin.itp.R")
source("use_for_all/FDR-R-Code/adpt.cutz.R")
source("use_for_all/FDR-R-Code/adaptZ.func.R")

#Sample size
n=200 #400
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

SEED<-scan("seed.txt")
seed = SEED[1] #1-100
set.seed(seed)
#Generate Gaussian data
data <- mvrnorm(n,Mean,Sigma2)
X_plus <- as.matrix(data)
X_plus2 <- matrixMultiply(t(X_plus),X_plus)/n

#Estimate precision matrix
j<-0
sigma<-matrix()
sigma_c<-matrix(1,p-1,1)

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
#Omega_hat_S
MIO_F_sigma1<-matrix_generator1(sigma,sigma_c) #Minimum symmetrization
#Omega_hat_US
MIO_F_sigma3<-matrix_generator3(sigma,sigma_c) #No symmetrization
TT <- matrixMultiply(t(MIO_F_sigma3),X_plus2)
#T_hat
MIO_T<-MIO_F_sigma3+t(MIO_F_sigma3)-matrixMultiply(TT,MIO_F_sigma3) #Desparsified estimator

MIO_F_lowertri1 <- MIO_F_sigma1
MIO_F_lowertri1[upper.tri(MIO_F_lowertri1)] <- 0
MIO_F_lowertri1_nodiag <- matrix(0,nrow=p-1,ncol=p)
MIO_F_lowertri1_nodiag[!upper.tri(MIO_F_lowertri1_nodiag)] <- MIO_F_lowertri1[lower.tri(MIO_F_lowertri1)]

#T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S))
j <- 0
zvalue_unbias_min <- {}
while(j<p){
  j<-j+1
  Ajhat <- which(MIO_F_lowertri1[-j,j]!=0) 
  nj <- length(Ajhat)
  Ajhat_plus <- which(MIO_F_sigma3[,j]!=0)
  jAjhat <- which(Ajhat_plus==j)
  if(nj==0){next}
  for(i in 1:nj){
    iAjhat <- which(Ajhat_plus==(Ajhat[i]+1))
    theta_hat_ij <- MIO_F_sigma3[-j,j][Ajhat][i] 
    sigmahat_ij_lowertri <- sqrt((solve(X_plus2[Ajhat_plus,Ajhat_plus])[iAjhat,iAjhat])*(solve(X_plus2[Ajhat_plus,Ajhat_plus])[jAjhat,jAjhat])+(solve(X_plus2[Ajhat_plus,Ajhat_plus])[iAjhat,jAjhat])^2)
    zvalue_ij <- sqrt(n)*theta_hat_ij/sigmahat_ij_lowertri
    zvalue_unbias_min <- c(zvalue_unbias_min,zvalue_ij)
  }
}
Ahat_unbias_min <- which(MIO_F_lowertri1_nodiag!=0)
critical_unbias_min <- Ahat_unbias_min[adaptZ.func(zvalue_unbias_min, alpha)$ac]
sigma_c_star_unbias_min <- MIO_F_lowertri1_nodiag
sigma_c_star_unbias_min[critical_unbias_min] <- 0
sigma_c_star_unbias_min[upper.tri(sigma_c_star_unbias_min)] <- t(sigma_c_star_unbias_min)[!lower.tri(t(sigma_c_star_unbias_min))]
#T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S))
MIO_T_star_unbias_min <- matrix_generator3(sigma,sigma_c_star_unbias_min) 

#T(Omega_hat_S|Z0(T_hat),SL(Omega_hat_S)) and T(T_hat|Z0(T_hat),SL(Omega_hat_S))
zvalue_desparsified_min <- computeZValueDesparsifiedGaussian(MIO_F_lowertri1,MIO_T,MIO_F_sigma3,sigma,n,p)
Ahat_desparsified_min <- which(MIO_F_lowertri1_nodiag!=0)
critical_desparsified_min <- Ahat_desparsified_min[adaptZ.func(zvalue_desparsified_min, alpha)$ac]
sigma_c_star_desparsified_min <- MIO_F_lowertri1_nodiag
sigma_c_star_desparsified_min[critical_desparsified_min] <- 0
sigma_c_star_desparsified_min[upper.tri(sigma_c_star_desparsified_min)] <- t(sigma_c_star_desparsified_min)[!lower.tri(t(sigma_c_star_desparsified_min))]
#T(Omega_hat_S|Z0(T_hat),SL(Omega_hat_S))
MIO_T_star_desparsified_min_method1 <- matrix_generator3(sigma,sigma_c_star_desparsified_min) 
Ahat_star_desparsified_min <- which(MIO_T_star_desparsified_min_method1==0)
#T(T_hat|Z0(T_hat),SL(Omega_hat_S))
MIO_T_star_desparsified_min_method2 <- MIO_T
MIO_T_star_desparsified_min_method2[Ahat_star_desparsified_min] <- 0 

#T(T_hat|Z0(T_hat),SL(T_hat))
zvalue_low <- computeZValueLowGaussian(MIO_T,MIO_F_sigma3,n,p)
adaptZ <- adaptZ.func(zvalue_low, alpha)
pvalue_low_adjust <- adaptZ$lfdr
pmatrix <- matrix(0,nrow=p,ncol=p)
pmatrix[lower.tri(pmatrix)] <- pvalue_low_adjust
pmatrix[upper.tri(pmatrix)] <- t(pmatrix)[upper.tri(t(pmatrix))]
critical_method3 <- which(pmatrix>adaptZ$threshold)
#T(T_hat|Z0(T_hat),SL(T_hat))
MIO_T_star_method3 <- MIO_T
MIO_T_star_method3[critical_method3] <- 0

#Support recovery performance
eva_min <- evaluation(omega2,MIO_F_sigma1)
eva_star_unbias_min <- evaluation(omega2,MIO_T_star_unbias_min)
eva_star_desparsified_min <- evaluation(omega2,MIO_T_star_desparsified_min_method1)
eva_star_method3 <- evaluation(omega2,MIO_T_star_method3)

precision_min <- eva_min$precision
precision_star_unbias_min <- eva_star_unbias_min$precision
sensitivity_min <- eva_min$sensitivity
sensitivity_star_unbias_min <- eva_star_unbias_min$sensitivity
specificity_min <- eva_min$specificity
specificity_star_unbias_min <- eva_star_unbias_min$specificity
MCC_min <- eva_min$MCC
MCC_star_unbias_min <- eva_star_unbias_min$MCC
F1_min <- eva_min$F1
F1_star_unbias_min <- eva_star_unbias_min$F1
precision_star_desparsified_min <- eva_star_desparsified_min$precision
sensitivity_star_desparsified_min <- eva_star_desparsified_min$sensitivity
specificity_star_desparsified_min <- eva_star_desparsified_min$specificity
MCC_star_desparsified_min <- eva_star_desparsified_min$MCC
F1_star_desparsified_min <- eva_star_desparsified_min$F1
precision_star_method3<- eva_star_method3$precision
sensitivity_star_method3 <- eva_star_method3$sensitivity
specificity_star_method3 <- eva_star_method3$specificity
MCC_star_method3 <- eva_star_method3$MCC
F1_star_method3 <- eva_star_method3$F1

#Matrix norm loss
MIO_min_norm1<-base::norm(MIO_F_sigma1-omega2,"1")
MIO_min_norm2<-base::norm(MIO_F_sigma1-omega2,"2")
MIO_min_normF<-base::norm(MIO_F_sigma1-omega2,"F")
MIO_min_normmax<-base::norm(MIO_F_sigma1-omega2,"M")
MIO_norm1_T_star_unbias_min<-base::norm(MIO_T_star_unbias_min-omega2,"1")
MIO_norm2_T_star_unbias_min<-base::norm(MIO_T_star_unbias_min-omega2,"2")
MIO_normF_T_star_unbias_min<-base::norm(MIO_T_star_unbias_min-omega2,"F")
MIO_normmax_T_star_unbias_min<-base::norm(MIO_T_star_unbias_min-omega2,"M")
MIO_norm1_T_star_desparsified_min_method1<-base::norm(MIO_T_star_desparsified_min_method1-omega2,"1")
MIO_norm2_T_star_desparsified_min_method1<-base::norm(MIO_T_star_desparsified_min_method1-omega2,"2")
MIO_normF_T_star_desparsified_min_method1<-base::norm(MIO_T_star_desparsified_min_method1-omega2,"F")
MIO_normmax_T_star_desparsified_min_method1<-base::norm(MIO_T_star_desparsified_min_method1-omega2,"M")
MIO_norm1_T_star_desparsified_min_method2<-base::norm(MIO_T_star_desparsified_min_method2-omega2,"1")
MIO_norm2_T_star_desparsified_min_method2<-base::norm(MIO_T_star_desparsified_min_method2-omega2,"2")
MIO_normF_T_star_desparsified_min_method2<-base::norm(MIO_T_star_desparsified_min_method2-omega2,"F")
MIO_normmax_T_star_desparsified_min_method2<-base::norm(MIO_T_star_desparsified_min_method2-omega2,"M")
MIO_norm1_T_star_method3<-base::norm(MIO_T_star_method3-omega2,"1")
MIO_norm2_T_star_method3<-base::norm(MIO_T_star_method3-omega2,"2")
MIO_normF_T_star_method3<-base::norm(MIO_T_star_method3-omega2,"F")
MIO_normmax_T_star_method3<-base::norm(MIO_T_star_method3-omega2,"M")
norm1_T <- base::norm(MIO_T-omega2,"1")
norm2_T <- base::norm(MIO_T-omega2,"2")
normF_T <- base::norm(MIO_T-omega2,"F")
normmax_T <- base::norm(MIO_T-omega2,"M")

save(precision_min,precision_star_unbias_min,
     sensitivity_min,sensitivity_star_unbias_min,specificity_min,
     specificity_star_unbias_min,MCC_min,MCC_star_unbias_min,F1_min,F1_star_unbias_min,
     precision_star_desparsified_min,
     sensitivity_star_desparsified_min,specificity_star_desparsified_min,
     MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
     sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
     MIO_min_norm1,MIO_min_norm2,MIO_min_normF,MIO_min_normmax,
     norm1_T,norm2_T,normF_T,normmax_T,
     MIO_norm1_T_star_unbias_min,MIO_norm2_T_star_unbias_min,
     MIO_normF_T_star_unbias_min,MIO_normmax_T_star_unbias_min,
     MIO_norm1_T_star_desparsified_min_method1,
     MIO_norm2_T_star_desparsified_min_method1,
     MIO_normF_T_star_desparsified_min_method1,
     MIO_normmax_T_star_desparsified_min_method1,
     MIO_norm1_T_star_desparsified_min_method2,
     MIO_norm2_T_star_desparsified_min_method2,
     MIO_normF_T_star_desparsified_min_method2,
     MIO_normmax_T_star_desparsified_min_method2,
     MIO_norm1_T_star_method3,MIO_norm2_T_star_method3,
     MIO_normF_T_star_method3,MIO_normmax_T_star_method3,file=paste("MIO_desparsified_unbias_Gaussian_band_n",n,"_p",p,"_seed",seed,"_one-two-stage_lfdr.RData",sep=""))

# save(precision_min,precision_star_unbias_min,
#      sensitivity_min,sensitivity_star_unbias_min,specificity_min,
#      specificity_star_unbias_min,MCC_min,MCC_star_unbias_min,F1_min,F1_star_unbias_min,
#      precision_star_desparsified_min,
#      sensitivity_star_desparsified_min,specificity_star_desparsified_min,
#      MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
#      sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
#      MIO_min_norm1,MIO_min_norm2,MIO_min_normF,MIO_min_normmax,
#      norm1_T,norm2_T,normF_T,normmax_T,
#      MIO_norm1_T_star_unbias_min,MIO_norm2_T_star_unbias_min,
#      MIO_normF_T_star_unbias_min,MIO_normmax_T_star_unbias_min,
#      MIO_norm1_T_star_desparsified_min_method1,
#      MIO_norm2_T_star_desparsified_min_method1,
#      MIO_normF_T_star_desparsified_min_method1,
#      MIO_normmax_T_star_desparsified_min_method1,
#      MIO_norm1_T_star_desparsified_min_method2,
#      MIO_norm2_T_star_desparsified_min_method2,
#      MIO_normF_T_star_desparsified_min_method2,
#      MIO_normmax_T_star_desparsified_min_method2,
#      MIO_norm1_T_star_method3,MIO_norm2_T_star_method3,
#      MIO_normF_T_star_method3,MIO_normmax_T_star_method3,file=paste("MIO_desparsified_unbias_Gaussian_random_n",n,"_p",p,"_seed",seed,"_one-two-stage_lfdr.RData",sep=""))
# 
# save(precision_min,precision_star_unbias_min,
#      sensitivity_min,sensitivity_star_unbias_min,specificity_min,
#      specificity_star_unbias_min,MCC_min,MCC_star_unbias_min,F1_min,F1_star_unbias_min,
#      precision_star_desparsified_min,
#      sensitivity_star_desparsified_min,specificity_star_desparsified_min,
#      MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
#      sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
#      MIO_min_norm1,MIO_min_norm2,MIO_min_normF,MIO_min_normmax,
#      norm1_T,norm2_T,normF_T,normmax_T,
#      MIO_norm1_T_star_unbias_min,MIO_norm2_T_star_unbias_min,
#      MIO_normF_T_star_unbias_min,MIO_normmax_T_star_unbias_min,
#      MIO_norm1_T_star_desparsified_min_method1,
#      MIO_norm2_T_star_desparsified_min_method1,
#      MIO_normF_T_star_desparsified_min_method1,
#      MIO_normmax_T_star_desparsified_min_method1,
#      MIO_norm1_T_star_desparsified_min_method2,
#      MIO_norm2_T_star_desparsified_min_method2,
#      MIO_normF_T_star_desparsified_min_method2,
#      MIO_normmax_T_star_desparsified_min_method2,
#      MIO_norm1_T_star_method3,MIO_norm2_T_star_method3,
#      MIO_normF_T_star_method3,MIO_normmax_T_star_method3,file=paste("MIO_desparsified_unbias_Gaussian_hub_n",n,"_p",p,"_seed",seed,"_one-two-stage_lfdr.RData",sep=""))
# 
# save(precision_min,precision_star_unbias_min,
#      sensitivity_min,sensitivity_star_unbias_min,specificity_min,
#      specificity_star_unbias_min,MCC_min,MCC_star_unbias_min,F1_min,F1_star_unbias_min,
#      precision_star_desparsified_min,
#      sensitivity_star_desparsified_min,specificity_star_desparsified_min,
#      MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
#      sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
#      MIO_min_norm1,MIO_min_norm2,MIO_min_normF,MIO_min_normmax,
#      norm1_T,norm2_T,normF_T,normmax_T,
#      MIO_norm1_T_star_unbias_min,MIO_norm2_T_star_unbias_min,
#      MIO_normF_T_star_unbias_min,MIO_normmax_T_star_unbias_min,
#      MIO_norm1_T_star_desparsified_min_method1,
#      MIO_norm2_T_star_desparsified_min_method1,
#      MIO_normF_T_star_desparsified_min_method1,
#      MIO_normmax_T_star_desparsified_min_method1,
#      MIO_norm1_T_star_desparsified_min_method2,
#      MIO_norm2_T_star_desparsified_min_method2,
#      MIO_normF_T_star_desparsified_min_method2,
#      MIO_normmax_T_star_desparsified_min_method2,
#      MIO_norm1_T_star_method3,MIO_norm2_T_star_method3,
#      MIO_normF_T_star_method3,MIO_normmax_T_star_method3,file=paste("MIO_desparsified_unbias_Gaussian_cluster_n",n,"_p",p,"_seed",seed,"_one-two-stage_lfdr.RData",sep=""))





