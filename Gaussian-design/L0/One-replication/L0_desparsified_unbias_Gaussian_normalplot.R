library(Rcpp)
library(RcppArmadillo)
library(huge)
library(MASS)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("use_for_all/asdar.cpp")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
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

#Estimate precision matrix
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
TT <- matrixMultiply(t(sdar_F_sigma),X_plus2)
#T_hat
sdar_T<-sdar_F_sigma+t(sdar_F_sigma)-matrixMultiply(TT,sdar_F_sigma) #Desparsified estimator

#Z-score of unbias estimator
z_unbias <- {}
j <- 0
while(j<p){
  j<-j+1
  Ajhat <- which(sdar_F_sigma[,j]!=0)
  nj <- length(Ajhat)
  jAjhat <- which(Ajhat==j)
  if(nj==0){next}
  for(i in 1:nj){
    theta_hat_ij <- sdar_F_sigma[,j][Ajhat][i]
    sigmahat_ij <- sqrt((solve(X_plus2[Ajhat,Ajhat])[i,i])*(solve(X_plus2[Ajhat,Ajhat])[jAjhat,jAjhat])+(solve(X_plus2[Ajhat,Ajhat])[i,jAjhat])^2)
    z_ij <- sqrt(n)*(theta_hat_ij-omega1[,j][Ajhat][i])/sigmahat_ij
    z_unbias <- c(z_unbias,z_ij)
  }
}

#Z-score of desparsified estimator
z_desparsified <- {}
j <- 0
while(j<p){
  j<-j+1
  Ajhat <- which(sdar_F_sigma[,j]!=0)
  nj <- length(Ajhat)
  if(nj==0){next}
  for(i in 1:nj){
    T_hat_ij <- sdar_T[,j][Ajhat][i]
    theta_hat_ij <- sdar_F_sigma[,j][Ajhat][i] 
    theta_hat_ji <- t(sdar_F_sigma)[,j][Ajhat][i] 
    sigmahat_ij2 <- as.matrix(sdar_F_sigma[Ajhat,Ajhat])[i,i]*sigma[j]+0.5*(theta_hat_ij^2+theta_hat_ji^2)
    sigmahat_ij <- sqrt(sigmahat_ij2)
    z_ij <- sqrt(n)*(T_hat_ij-omega1[,j][Ajhat][i])/sigmahat_ij
    z_desparsified <- c(z_desparsified,z_ij)
  }
}

save(Ahat,z_unbias,z_desparsified,file=paste("L0_desparsified_unbias_Gaussian_band_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))

# save(Ahat,z_unbias,z_desparsified,file=paste("L0_desparsified_unbias_Gaussian_random_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))
# 
# save(Ahat,z_unbias,z_desparsified,file=paste("L0_desparsified_unbias_Gaussian_hub_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))
# 
# save(Ahat,z_unbias,z_desparsified,file=paste("L0_desparsified_unbias_Gaussian_cluster_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))


