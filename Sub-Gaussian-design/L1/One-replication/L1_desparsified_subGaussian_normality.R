library(Rcpp)
library(RcppArmadillo)
library(huge)
library(MASS)
library(glmnet)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("use_for_all/asdar.cpp")
source("use_for_all/ans-end.R")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
source("use_for_all/useful_function.R")

#Sample size
n=200 #400
#Dimension
p=200 #400,1000,4000
#FDR level
alpha <- 0.05

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

#Candidate values for tunning parameter
num_lambda=20;lambda_max=2;b=50;lambda_min_ratio=0.01
lambdamax = lambda_max
nlambda=num_lambda
lambda.min.ratio=lambda_min_ratio
lambda.min = lambda.min.ratio*lambdamax
Lambda = exp(seq(log(lambdamax), log(lambda.min), length = nlambda))

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

result <- normalizeMatrix(data)
data_normalize <- result$normalizedX
diag <- result$norms

#Estimate support fot L0
j<-0
Beta <- {}
while(j<p){
  j<-j+1
  X<-data[,-j]
  y<-data[,j]
  nX <- data_normalize[,-j]
  dx <- diag[-j]
  # X<-as.matrix(X) 
  y<-as.vector(y)
  best_ebeta<-asdar(X,nX,dx,y,20,n,p)
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
TT <- matrixMultiply(t(NS_F_sigma),X_plus2)
#T_hat
NS_T<-NS_F_sigma+t(NS_F_sigma)-matrixMultiply(TT,NS_F_sigma) #Desparsified estimator

dataSigma <- matrixMultiply(data,NS_F_sigma)
dataSigma2 <- dataSigma^2

result <- L1_desparsified_confidence_intervals_SubGaussian(NS_T,NS_F_sigma,omega2,dataSigma2,n,p,alpha)
cov <- result$cov_desparsified #Coverage rate
z <- result$z_desparsified #Z-score
l <- result$l_desparsified #Avglength

z <- round(z,3)
l <- round(l,3)

save(Ahat,cov,z,l,file=paste("L1_desparsified_subGaussian_band_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov,z,l,file=paste("L1_desparsified_subGaussian_random_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov,z,l,file=paste("L1_desparsified_subGaussian_hub_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov,z,l,file=paste("L1_desparsified_subGaussian_cluster_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))