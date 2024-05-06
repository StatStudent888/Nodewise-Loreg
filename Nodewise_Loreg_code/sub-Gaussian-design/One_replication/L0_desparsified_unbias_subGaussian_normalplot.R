library(Rcpp)
library(RcppArmadillo)
library(huge)
library(MASS)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("end-sdar/sdar_RCPP-end.cpp")
source("use_for_all/asdar-end.R")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")

#Sample size
n=200 #400,800
#Dimension
p=200 #400

#Generate Band graph
set.seed(-1)
omega1 <- omega_generator1(p,1,0.5,0.3)
Sigma1 <- solve(omega1)

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
seed = SEED[1]

Sigma_sqrt <- eigen(Sigma1)$vectors%*%sqrt(diag(eigen(Sigma1)$values))%*%solve(eigen(Sigma1)$vectors)
U <- matrix(0,n,p)
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
Beta <- {}
sigma<-matrix()
sigma_c<-matrix(1,p-1,1)
while(j<p){
  j<-j+1
  X<-data[,-j]
  y<-data[,j]
  X<-as.matrix(X)
  y<-as.vector(y)
  best_ebeta<-asdar(X,y,20,n,p)
  Beta<-cbind(Beta,best_ebeta)
  sigmak<-n/(norm(y-X%*%best_ebeta,"2")^2)  
  sigma_ck<-(-1)*sigmak*best_ebeta   
  sigma<- cbind(sigma,sigmak)
  sigma_c<- cbind(sigma_c,sigma_ck)
}
Ahat <- which(Beta!=0)
sigma<-sigma[-1]
sigma_c<-sigma_c[,-1]
#Omega_hat_US
sdar_F_sigma<-matrix_generator3(sigma,sigma_c) #No symmetrization
#T_hat
sdar_T<-sdar_F_sigma+t(sdar_F_sigma)-t(sdar_F_sigma)%*%(t(data)%*%data/n)%*%sdar_F_sigma #Desparsified estimator

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
    theta_hat_ji <- t(sdar_F_sigma)[,j][Ajhat][i]
    B <- solve((t(X_plus)%*%X_plus/n)[Ajhat,Ajhat])
    sigmahat_ij <- sqrt((((B[i,])%*%t(as.matrix(X_plus[,Ajhat])))^2%*%((as.matrix(X_plus[,Ajhat]))%*%(as.matrix(B[,jAjhat])))^2)/n-(1/2)*(theta_hat_ij^2+theta_hat_ji^2))
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
    theta_hat_ij <- sdar_F_sigma[,j][Ajhat][i] #代入未对称化的theta_hat_ij
    theta_hat_ji <- t(sdar_F_sigma)[,j][Ajhat][i] #代入未对称化的theta_hat_ji
    sigmahat_ij <- sqrt(((t(as.matrix(as.matrix(sdar_F_sigma[,Ajhat])[,i]))%*%t(data))^2%*%(data%*%as.matrix(sdar_F_sigma[,j]))^2)/n-(1/2)*(theta_hat_ij^2+theta_hat_ji^2))
    z_ij <- sqrt(n)*(T_hat_ij-omega1[,j][Ajhat][i])/sigmahat_ij
    z_desparsified <- c(z_desparsified,z_ij)
  }
}

save(Ahat,z_unbias,z_desparsified,sdar_F_sigma,sdar_T,data,file=paste("L0_desparsified_unbias_subGaussian_band_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))

#save(Ahat,z_unbias,z_desparsified,sdar_F_sigma,sdar_T,data,file=paste("L0_desparsified_unbias_subGaussian_random_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))

#save(Ahat,z_unbias,z_desparsified,sdar_F_sigma,sdar_T,data,file=paste("L0_desparsified_unbias_subGaussian_hub_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))

#save(Ahat,z_unbias,z_desparsified,sdar_F_sigma,sdar_T,data,file=paste("L0_desparsified_unbias_subGaussian_cluster_n",n,"_p",p,"_seed",seed,"_normalplot.RData",sep=""))






