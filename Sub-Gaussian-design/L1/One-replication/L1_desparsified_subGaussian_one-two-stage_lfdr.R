library(Rcpp)
library(RcppArmadillo)
library(huge)
library(glmnet)
sourceCpp("use_for_all/asdar.cpp")
source("use_for_all/ans-end.R")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
source("use_for_all/evaluation.R")
source("use_for_all/FDR-R-Code/epsest.func.R")
source("use_for_all/FDR-R-Code/lin.itp.R")
source("use_for_all/FDR-R-Code/adpt.cutz.R")
source("use_for_all/FDR-R-Code/adaptZ.func.R")

#Sample size
n=200 #400
#Dimension
p=200 #400,1000,4000
#FDR level
alpha <- 0.05

#Candidate values for tunning parameter
num_lambda=20;lambda_max=2;b=50;lambda_min_ratio=0.01
lambdamax = lambda_max
nlambda=num_lambda
lambda.min.ratio=lambda_min_ratio
lambda.min = lambda.min.ratio*lambdamax
Lambda = exp(seq(log(lambdamax), log(lambda.min), length = nlambda))

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
  best_ebeta_lambda<-ans(X,y,Lambda,n,p)
  best_ebeta <- best_ebeta_lambda$bestbeta
  best_lambda <- best_ebeta_lambda$bestlambda
  sigmak<-1/(((norm(y-X%*%best_ebeta,"2")^2)/n)+best_lambda/2*sum(abs(best_ebeta))) 
  sigma_ck<-(-1)*sigmak*best_ebeta  
  sigma<- cbind(sigma,sigmak)
  sigma_c<- cbind(sigma_c,sigma_ck)
}
sigma<-sigma[-1]
sigma_c<-sigma_c[,-1]
#Omega_hat_S
NS_F_sigma1<-matrix_generator1(sigma,sigma_c) #Minimum symmetrization
#Omega_hat_US
NS_F_sigma3<-matrix_generator3(sigma,sigma_c) #No symmetrization
TT <- matrixMultiply(t(NS_F_sigma3),X_plus2)
#T_hat
NS_T<-NS_F_sigma3+t(NS_F_sigma3)-matrixMultiply(TT,NS_F_sigma3) #Desparsified estimator

NS_F_lowertri1 <- NS_F_sigma1
NS_F_lowertri1[upper.tri(NS_F_lowertri1)] <- 0
NS_F_lowertri1_nodiag <- matrix(0,nrow=p-1,ncol=p)
NS_F_lowertri1_nodiag[!upper.tri(NS_F_lowertri1_nodiag)] <- NS_F_lowertri1[lower.tri(NS_F_lowertri1)]

#T(Omega_hat_S|Z0(T_hat),SL(Omega_hat_S)) and T(T_hat|Z0(T_hat),SL(Omega_hat_S))
dataSigma <- matrixMultiply(data,NS_F_sigma3)
dataSigma2 <- dataSigma^2
j <- 0
zvalue_desparsified_min <- {}
while(j<p){
  j<-j+1
  Ajhat <- which(NS_F_lowertri1[-j,j]!=0) 
  nj <- length(Ajhat)
  if(nj==0){next}
  for(i in 1:nj){
    T_hat_ij <- NS_T[-j,j][Ajhat][i]
    theta_hat_ij <- NS_F_sigma3[-j,j][Ajhat][i] 
    theta_hat_ji <- t(NS_F_sigma3)[-j,j][Ajhat][i] 
    sigmahat_ij2 <- t(as.matrix(as.matrix(dataSigma2[,-j][,Ajhat])[,i]))%*%as.matrix(dataSigma2[,j])/n-(0.5*(theta_hat_ij^2+theta_hat_ji^2))
    sigmahat_ij_lowertri <- sqrt(sigmahat_ij2)
    zvalue_ij <- sqrt(n)*T_hat_ij/sigmahat_ij_lowertri
    zvalue_desparsified_min <- c(zvalue_desparsified_min,zvalue_ij)
  }
}
Ahat_desparsified_min <- which(NS_F_lowertri1_nodiag!=0)
critical_desparsified_min <- Ahat_desparsified_min[adaptZ.func(zvalue_desparsified_min, alpha)$ac]
sigma_c_star_desparsified_min <- NS_F_lowertri1_nodiag
sigma_c_star_desparsified_min[critical_desparsified_min] <- 0
sigma_c_star_desparsified_min[upper.tri(sigma_c_star_desparsified_min)] <- t(sigma_c_star_desparsified_min)[!lower.tri(t(sigma_c_star_desparsified_min))]
#T(Omega_hat_S|Z0(T_hat),SL(Omega_hat_S))
NS_T_star_desparsified_min_method1 <- matrix_generator3(sigma,sigma_c_star_desparsified_min) 
Ahat_star_desparsified_min <- which(NS_T_star_desparsified_min_method1==0)
#T(T_hat|Z0(T_hat),SL(Omega_hat_S))
NS_T_star_desparsified_min_method2 <- NS_T
NS_T_star_desparsified_min_method2[Ahat_star_desparsified_min] <- 0 

#T(T_hat|Z0(T_hat),SL(T_hat))
zvalue_low <- computeZValueLowSubGaussian(NS_T,NS_F_sigma3,dataSigma2,n,p)
adaptZ <- adaptZ.func(zvalue_low, alpha)
pvalue_low_adjust <- adaptZ$lfdr
pmatrix <- matrix(0,nrow=p,ncol=p)
pmatrix[lower.tri(pmatrix)] <- pvalue_low_adjust
pmatrix[upper.tri(pmatrix)] <- t(pmatrix)[upper.tri(t(pmatrix))]
critical_method3 <- which(pmatrix>adaptZ$threshold)
#T(T_hat|Z0(T_hat),SL(T_hat))
NS_T_star_method3 <- NS_T
NS_T_star_method3[critical_method3] <- 0

#Support recovery performance
eva_min <- evaluation(omega2,NS_F_sigma1)
eva_star_desparsified_min <- evaluation(omega2,NS_T_star_desparsified_min_method1)
eva_star_method3 <- evaluation(omega2,NS_T_star_method3)

precision_min <- eva_min$precision
sensitivity_min <- eva_min$sensitivity
specificity_min <- eva_min$specificity
MCC_min <- eva_min$MCC
F1_min <- eva_min$F1
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
NS_min_norm1<-norm(NS_F_sigma1-omega2,"1")
NS_min_norm2<-norm(NS_F_sigma1-omega2,"2")
NS_min_normF<-norm(NS_F_sigma1-omega2,"F")
NS_min_normmax<-norm(NS_F_sigma1-omega2,"M")
NS_norm1_T_star_desparsified_min_method1<-norm(NS_T_star_desparsified_min_method1-omega2,"1")
NS_norm2_T_star_desparsified_min_method1<-norm(NS_T_star_desparsified_min_method1-omega2,"2")
NS_normF_T_star_desparsified_min_method1<-norm(NS_T_star_desparsified_min_method1-omega2,"F")
NS_normmax_T_star_desparsified_min_method1<-norm(NS_T_star_desparsified_min_method1-omega2,"M")
NS_norm1_T_star_desparsified_min_method2<-norm(NS_T_star_desparsified_min_method2-omega2,"1")
NS_norm2_T_star_desparsified_min_method2<-norm(NS_T_star_desparsified_min_method2-omega2,"2")
NS_normF_T_star_desparsified_min_method2<-norm(NS_T_star_desparsified_min_method2-omega2,"F")
NS_normmax_T_star_desparsified_min_method2<-norm(NS_T_star_desparsified_min_method2-omega2,"M")
NS_norm1_T_star_method3<-norm(NS_T_star_method3-omega2,"1")
NS_norm2_T_star_method3<-norm(NS_T_star_method3-omega2,"2")
NS_normF_T_star_method3<-norm(NS_T_star_method3-omega2,"F")
NS_normmax_T_star_method3<-norm(NS_T_star_method3-omega2,"M")
norm1_T <- norm(NS_T-omega2,"1")
norm2_T <- norm(NS_T-omega2,"2")
normF_T <- norm(NS_T-omega2,"F")
normmax_T <- norm(NS_T-omega2,"M")

save(precision_min,
     sensitivity_min,specificity_min,
     MCC_min,F1_min,
     precision_star_desparsified_min,
     sensitivity_star_desparsified_min,specificity_star_desparsified_min,
     MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
     sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
     NS_min_norm1,NS_min_norm2,NS_min_normF,NS_min_normmax,
     norm1_T,norm2_T,normF_T,normmax_T,
     NS_norm1_T_star_desparsified_min_method1,
     NS_norm2_T_star_desparsified_min_method1,
     NS_normF_T_star_desparsified_min_method1,
     NS_normmax_T_star_desparsified_min_method1,
     NS_norm1_T_star_desparsified_min_method2,
     NS_norm2_T_star_desparsified_min_method2,
     NS_normF_T_star_desparsified_min_method2,
     NS_normmax_T_star_desparsified_min_method2,
     NS_norm1_T_star_method3,NS_norm2_T_star_method3,
     NS_normF_T_star_method3,NS_normmax_T_star_method3,file=paste("L1_desparsified_subGaussian_band_n",n,"_p",p,"_seed",seed,"_one-two-stage_lfdr.RData",sep=""))

# save(precision_min,
#      sensitivity_min,specificity_min,
#      MCC_min,F1_min,
#      precision_star_desparsified_min,
#      sensitivity_star_desparsified_min,specificity_star_desparsified_min,
#      MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
#      sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
#      NS_min_norm1,NS_min_norm2,NS_min_normF,NS_min_normmax,
#      norm1_T,norm2_T,normF_T,normmax_T,
#      NS_norm1_T_star_desparsified_min_method1,
#      NS_norm2_T_star_desparsified_min_method1,
#      NS_normF_T_star_desparsified_min_method1,
#      NS_normmax_T_star_desparsified_min_method1,
#      NS_norm1_T_star_desparsified_min_method2,
#      NS_norm2_T_star_desparsified_min_method2,
#      NS_normF_T_star_desparsified_min_method2,
#      NS_normmax_T_star_desparsified_min_method2,
#      NS_norm1_T_star_method3,NS_norm2_T_star_method3,
#      NS_normF_T_star_method3,NS_normmax_T_star_method3,file=paste("L1_desparsified_subGaussian_random_n",n,"_p",p,"_seed",seed,"_one-two-stage_lfdr.RData",sep=""))
# 
# save(precision_min,
#      sensitivity_min,specificity_min,
#      MCC_min,F1_min,
#      precision_star_desparsified_min,
#      sensitivity_star_desparsified_min,specificity_star_desparsified_min,
#      MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
#      sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
#      NS_min_norm1,NS_min_norm2,NS_min_normF,NS_min_normmax,
#      norm1_T,norm2_T,normF_T,normmax_T,
#      NS_norm1_T_star_desparsified_min_method1,
#      NS_norm2_T_star_desparsified_min_method1,
#      NS_normF_T_star_desparsified_min_method1,
#      NS_normmax_T_star_desparsified_min_method1,
#      NS_norm1_T_star_desparsified_min_method2,
#      NS_norm2_T_star_desparsified_min_method2,
#      NS_normF_T_star_desparsified_min_method2,
#      NS_normmax_T_star_desparsified_min_method2,
#      NS_norm1_T_star_method3,NS_norm2_T_star_method3,
#      NS_normF_T_star_method3,NS_normmax_T_star_method3,file=paste("L1_desparsified_subGaussian_hub_n",n,"_p",p,"_seed",seed,"_one-two-stage_lfdr.RData",sep=""))
# 
# save(precision_min,
#      sensitivity_min,specificity_min,
#      MCC_min,F1_min,
#      precision_star_desparsified_min,
#      sensitivity_star_desparsified_min,specificity_star_desparsified_min,
#      MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
#      sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
#      NS_min_norm1,NS_min_norm2,NS_min_normF,NS_min_normmax,
#      norm1_T,norm2_T,normF_T,normmax_T,
#      NS_norm1_T_star_desparsified_min_method1,
#      NS_norm2_T_star_desparsified_min_method1,
#      NS_normF_T_star_desparsified_min_method1,
#      NS_normmax_T_star_desparsified_min_method1,
#      NS_norm1_T_star_desparsified_min_method2,
#      NS_norm2_T_star_desparsified_min_method2,
#      NS_normF_T_star_desparsified_min_method2,
#      NS_normmax_T_star_desparsified_min_method2,
#      NS_norm1_T_star_method3,NS_norm2_T_star_method3,
#      NS_normF_T_star_method3,NS_normmax_T_star_method3,file=paste("L1_desparsified_subGaussian_cluster_n",n,"_p",p,"_seed",seed,"_one-two-stage_lfdr.RData",sep=""))








