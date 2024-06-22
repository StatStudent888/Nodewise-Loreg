library(Rcpp)
library(RcppArmadillo)
library(huge)
library(MASS)
library(glmnet)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("end-sdar/sdar_RCPP-end.cpp")
source("use_for_all/asdar-end.R")
source("use_for_all/ans-end.R")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")

#Sample size
n=200 #400
#Dimension
p=200 #400
#FDR level
alpha <- 0.05

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

#Candidate values for tunning parameter
num_lambda=20;lambda_max=2;b=50;lambda_min_ratio=0.01
lambdamax = lambda_max
nlambda=num_lambda
lambda.min.ratio=lambda_min_ratio
lambda.min = lambda.min.ratio*lambdamax
Lambda = exp(seq(log(lambdamax), log(lambda.min), length = nlambda))

SEED<-scan("seed.txt")
seed = SEED[1]
set.seed(seed)
#Generate Gaussian data
data <- mvrnorm(n,Mean,Sigma2)

#Get the support of Nodewise Loreg algorithm
j<-0
Beta <- {}
while(j<p){
  j<-j+1
  X<-data[,-j]
  y<-data[,j]
  X<-as.matrix(X) 
  y<-as.vector(y)
  best_ebeta<-asdar(X,y,20,n,p)
  Beta<-cbind(Beta,best_ebeta)
}
Ahat <- which(Beta!=0)

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
#T_hat
NS_T<-NS_F_sigma+t(NS_F_sigma)-t(NS_F_sigma)%*%(t(data)%*%data/n)%*%NS_F_sigma #Desparsified estimator

cov <- {} #Coverage rate
z <- {} #Z-score
l <- {} #Avglength
for(j in 1:p){
  for(i in 1:p){
    T_hat_ij <- NS_T[i,j]
    sigmahat_ij2 <- NS_F_sigma[i,i]*NS_F_sigma[j,j]+(NS_F_sigma[i,j])^2
    sigmahat_ji2 <- NS_F_sigma[i,i]*NS_F_sigma[j,j]+(NS_F_sigma[j,i])^2
    sigmahat_ij <- sqrt((sigmahat_ij2+sigmahat_ji2)/2)
    z_ij <- sqrt(n)*(T_hat_ij-omega2[i,j])/sigmahat_ij
    z <- c(z,z_ij)
    l_ij <- 2*qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
    l <- c(l,l_ij)
    CI_ij_low <- T_hat_ij-qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
    CI_ij_up <- T_hat_ij+qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
    if((omega2[i,j])<=CI_ij_up & (omega2[i,j])>=CI_ij_low){
      cov_ij <- 1
    }
    else{cov_ij <- 0}
    cov <- c(cov,cov_ij)
  }
}

save(Ahat,cov,z,l,file=paste("L1_desparsified_Gaussian_band_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov,z,l,file=paste("L1_desparsified_Gaussian_random_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov,z,l,file=paste("L1_desparsified_Gaussian_hub_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov,z,l,file=paste("L1_desparsified_Gaussian_cluster_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))