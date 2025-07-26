library(Rcpp)
library(RcppArmadillo)
library(huge)
library(MASS)
library(glmnet)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("use_for_all/asdar.cpp")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
source("use_for_all/ans-end.R")
source("use_for_all/useful_function.R")

#Sample size
n=200 #400,800
#Dimension
p=200 #400,1000,4000

#Generate Band graph
set.seed(-1)
omega1 <- omega_generator1(p,1,0.5,0.3)
if(p %in% c(200,400)){
  Sigma1 <- solve(omega1)
} else {Sigma1 <- matrix_inverse(omega1)}

#Generate Random graph
#graph= "random"
#set.seed(-1)
#generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, prob=4/p)
#omega1 <- generator$omega
#Sigma1 <- generator$sigma

#Generate Hub graph
#graph= "hub"
#set.seed(-1)
#generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
#omega1 <- generator$omega
#Sigma1 <- generator$sigma

#Generate Cluster graph
#graph= "cluster"
#set.seed(-1)
#generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
#omega1 <- generator$omega
#Sigma1 <- generator$sigma

Mean <- rep(0,p)

#Candidate values for tunning parameter
num_lambda=20;lambda_max=2;b=50;lambda_min_ratio=0.01
lambdamax = lambda_max
nlambda=num_lambda
lambda.min.ratio=lambda_min_ratio
lambda.min = lambda.min.ratio*lambdamax
Lambda = exp(seq(log(lambdamax), log(lambda.min), length = nlambda))

SEED<-1:400
seed = SEED[1] #1-400
set.seed(seed)
#Generate Gaussian data
data <- mvrnorm(n,Mean,Sigma1)
X_plus <- as.matrix(data)
X_plus2 <- matrixMultiply(t(X_plus),X_plus)/n

result <- normalizeMatrix(data)
data_normalize <- result$normalizedX
diag <- result$norms

#Estimate support for L0
j<-0
Beta <- {}
sigma<-matrix()
sigma_c<-matrix(1,p-1,1)
while(j<p){
  j<-j+1
  X<-data[,-j]
  y<-data[,j]
  nX <- data_normalize[,-j]
  dx <- diag[-j]
  # X<-as.matrix(X) 
  y<-as.vector(y)
  best_ebeta<-asdar(X,nX,dx,y,20,n,p)
  sigmak<-n/(norm(y-X%*%best_ebeta,"2")^2)  
  sigma_ck<-(-1)*sigmak*best_ebeta   
  sigma<- cbind(sigma,sigmak)
  sigma_c<- cbind(sigma_c,sigma_ck)
  Beta<-cbind(Beta,best_ebeta)
}
Ahat <- which(Beta!=0)
sigma<-sigma[-1]
sigma_c<-sigma_c[,-1]
#Omega_hat_US
sdar_F_sigma<-matrix_generator3(sigma,sigma_c) #No symmetrization

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
#Omega_hat_US
NS_F_sigma<-matrix_generator3(sigma,sigma_c) #No symmetrization
TT <- matrixMultiply(t(NS_F_sigma),X_plus2)
#T_hat
NS_T<-NS_F_sigma+t(NS_F_sigma)-matrixMultiply(TT,NS_F_sigma) #Desparsified estimator

#Z-score of desparsified estimator
z <- {}
j <- 0
while(j<p){
  j<-j+1
  Ajhat <- which(sdar_F_sigma[,j]!=0)
  nj <- length(Ajhat)
  if(nj==0){next}
  for(i in 1:nj){
    T_hat_ij <- NS_T[,j][Ajhat][i]
    theta_hat_ij <- NS_F_sigma[,j][Ajhat][i] 
    theta_hat_ji <- t(NS_F_sigma)[,j][Ajhat][i] 
    sigmahat_ij2 <- as.matrix(NS_F_sigma[Ajhat,Ajhat])[i,i]*sigma[j]+0.5*(theta_hat_ij^2+theta_hat_ji^2)
    sigmahat_ij <- sqrt(sigmahat_ij2)
    z_ij <- sqrt(n)*(T_hat_ij-omega1[,j][Ajhat][i])/sigmahat_ij
    z <- c(z,z_ij)
  }
}

save(Ahat,z,file=paste("L1_desparsified_Gaussian_band_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))

# save(Ahat,z,file=paste("L1_desparsified_Gaussian_random_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))
# 
# save(Ahat,z,file=paste("L1_desparsified_Gaussian_hub_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))
# 
# save(Ahat,z,file=paste("L1_desparsified_Gaussian_cluster_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))



