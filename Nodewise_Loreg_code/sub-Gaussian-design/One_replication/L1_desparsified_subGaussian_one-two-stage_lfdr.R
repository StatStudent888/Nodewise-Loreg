library(huge)
library(glmnet)
source("use_for_all/ans-end.R")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
source("use_for_all/evaluation.R")
source("FDR-R-Code/epsest.func.R")
source("FDR-R-Code/lin.itp.R")
source("FDR-R-Code/adpt.cutz.R")
source("FDR-R-Code/adaptZ.func.R")

#Sample size
n=200 #400
#Dimension
p=200 #400
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
Sigma2 <- solve(omega2)

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

Sigma_sqrt <- eigen(Sigma2)$vectors%*%sqrt(diag(eigen(Sigma2)$values))%*%solve(eigen(Sigma2)$vectors)
U <- matrix(0,n,p)

SEED<-scan("seed.txt")
seed = SEED[1]
set.seed(seed)

for(h in 1:n){
  for(q in 1:p){
    U[h,q] <- runif(1,-sqrt(3),sqrt(3))
  }
}
#Generate sub-Gaussian data
data <- U%*%Sigma_sqrt
X_plus <- as.matrix(data)

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
NS_F_sigma1<-matrix_generator1(sigma,sigma_c) 
#Omega_hat_US
NS_F_sigma3<-matrix_generator3(sigma,sigma_c) #No symmetrization
#T_hat
NS_T<-NS_F_sigma3+t(NS_F_sigma3)-t(NS_F_sigma3)%*%(t(data)%*%data/n)%*%NS_F_sigma3 #Desparsified estimator

NS_F_lowertri1 <- NS_F_sigma1
NS_F_lowertri1[upper.tri(NS_F_lowertri1)] <- 0
NS_F_lowertri1_nodiag <- matrix(0,nrow=p-1,ncol=p)
NS_F_lowertri1_nodiag[!upper.tri(NS_F_lowertri1_nodiag)] <- NS_F_lowertri1[lower.tri(NS_F_lowertri1)]

#T(Omega_hat_S|Z0(T_hat),SL(Omega_hat_S)) and T(T_hat|Z0(T_hat),SL(Omega_hat_S))
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
    sigmahat_ij2 <- ((t(as.matrix(as.matrix(NS_F_sigma3[,-j][,Ajhat])[,i]))%*%t(data))^2%*%(data%*%NS_F_sigma3[,j])^2)/n-theta_hat_ij^2
    sigmahat_ji2 <- ((t(as.matrix(NS_F_sigma3[,j]))%*%t(data))^2%*%(data%*%as.matrix(as.matrix(NS_F_sigma3[,-j][,Ajhat])[,i]))^2)/n-theta_hat_ji^2
    sigmahat_ij_lowertri <- sqrt((sigmahat_ij2+sigmahat_ji2)/2)
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
NS_T_star_desparsified_min_method1 <- matrix_generator3(sigma,sigma_c_star_desparsified_min) #绝对值最小对称化
Ahat_star_desparsified_min <- which(NS_T_star_desparsified_min_method1==0)
#T(T_hat|Z0(T_hat),SL(Omega_hat_S))
NS_T_star_desparsified_min_method2 <- NS_T
NS_T_star_desparsified_min_method2[Ahat_star_desparsified_min] <- 0 

#T(T_hat|Z0(T_hat),SL(T_hat))
zvalue_low <- {}
for(j in 1:p){
  for(i in 1:p){
    if(i>j){
      T_hat_ij <- NS_T[i,j]
      sigmahat_ij2 <- ((t(as.matrix(NS_F_sigma3[,i]))%*%t(data))^2%*%(data%*%NS_F_sigma3[,j])^2)/n-(NS_F_sigma3[i,j])^2
      sigmahat_ji2 <- ((t(as.matrix(NS_F_sigma3[,j]))%*%t(data))^2%*%(data%*%NS_F_sigma3[,i])^2)/n-(NS_F_sigma3[j,i])^2
      sigmahat_ij_lowertri <- sqrt((sigmahat_ij2+sigmahat_ji2)/2)
      zvalue_low_ij <- sqrt(n)*T_hat_ij/sigmahat_ij_lowertri
      zvalue_low <- c(zvalue_low,zvalue_low_ij)
    }
    else if(i<j | i==j){
      next
    }
  }
}
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

save(data,NS_F_sigma1,NS_F_sigma3,NS_T,
     NS_T_star_desparsified_min_method1,NS_T_star_desparsified_min_method2,
     NS_T_star_method3,precision_min,
     sensitivity_min,specificity_min,
     MCC_min,F1_min,
     precision_star_desparsified_min,
     sensitivity_star_desparsified_min,specificity_star_desparsified_min,
     MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
     sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
     NS_min_norm1,NS_min_norm2,NS_min_normF,NS_min_normmax,
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

# save(data,NS_F_sigma1,NS_F_sigma3,NS_T,
#      NS_T_star_desparsified_min_method1,NS_T_star_desparsified_min_method2,
#      NS_T_star_method3,precision_min,
#      sensitivity_min,specificity_min,
#      MCC_min,F1_min,
#      precision_star_desparsified_min,
#      sensitivity_star_desparsified_min,specificity_star_desparsified_min,
#      MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
#      sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
#      NS_min_norm1,NS_min_norm2,NS_min_normF,NS_min_normmax,
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

# save(data,NS_F_sigma1,NS_F_sigma3,NS_T,
#      NS_T_star_desparsified_min_method1,NS_T_star_desparsified_min_method2,
#      NS_T_star_method3,precision_min,
#      sensitivity_min,specificity_min,
#      MCC_min,F1_min,
#      precision_star_desparsified_min,
#      sensitivity_star_desparsified_min,specificity_star_desparsified_min,
#      MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
#      sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
#      NS_min_norm1,NS_min_norm2,NS_min_normF,NS_min_normmax,
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

# save(data,NS_F_sigma1,NS_F_sigma3,NS_T,
#      NS_T_star_desparsified_min_method1,NS_T_star_desparsified_min_method2,
#      NS_T_star_method3,precision_min,
#      sensitivity_min,specificity_min,
#      MCC_min,F1_min,
#      precision_star_desparsified_min,
#      sensitivity_star_desparsified_min,specificity_star_desparsified_min,
#      MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
#      sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
#      NS_min_norm1,NS_min_norm2,NS_min_normF,NS_min_normmax,
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