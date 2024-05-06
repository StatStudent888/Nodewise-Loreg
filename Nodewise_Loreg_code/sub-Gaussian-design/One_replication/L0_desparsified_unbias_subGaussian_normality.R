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
n=200 #400
#Dimension
p=200 #400

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
#T_hat
sdar_T<-sdar_F_sigma+t(sdar_F_sigma)-t(sdar_F_sigma)%*%(t(data)%*%data/n)%*%sdar_F_sigma #Desparsified estimator

cov_desparsified <- {} #Coverage rate
z_desparsified <- {} #Z-score
l_desparsified <- {} #Avglength
ltrue_desparsified <- {} #Truelength
for(j in 1:p){
  for(i in 1:p){
    T_hat_ij <- sdar_T[i,j]
    A <- Sigma_sqrt%*%((1/2)*(as.matrix(omega2[,i])%*%omega2[j,]+as.matrix(omega2[,j])%*%omega2[i,]))%*%Sigma_sqrt
    sigmahat_ij <- sqrt(((t(as.matrix(sdar_F_sigma[,i]))%*%t(data))^2%*%(data%*%as.matrix(sdar_F_sigma[,j]))^2)/n-(1/2)*((sdar_F_sigma[i,j])^2+(sdar_F_sigma[j,i])^2))
    sigmatrue_ij <- sqrt(((2*(sqrt(3))^4)/9)*(sum(diag(A%*%A)))-((2*(sqrt(3))^4)/15)*(sum((diag(A))^2)))
    z_ij <- sqrt(n)*(T_hat_ij-omega2[i,j])/sigmahat_ij
    z_desparsified <- c(z_desparsified,z_ij)
    l_ij <- 2*qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
    l_desparsified <- c(l_desparsified,l_ij)
    ltrue_ij <- 2*qnorm(1-alpha/2)*sigmatrue_ij/sqrt(n)
    ltrue_desparsified <- c(ltrue_desparsified,ltrue_ij)
    CI_ij_low <- T_hat_ij-qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
    CI_ij_up <- T_hat_ij+qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
    if((omega2[i,j])<=CI_ij_up & (omega2[i,j])>=CI_ij_low){
      cov_ij <- 1
    }
    else{cov_ij <- 0}
    cov_desparsified <- c(cov_desparsified,cov_ij)
  }
}

cov_unbias <- {} #Coverage rate
z_unbias <- {} #Z-score
l_unbias <- {} #Avglength
ltrue_unbias <- {} #Truelength
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
    B_true <- solve(Sigma2[Ajhat,Ajhat])
    A <- as.matrix(Sigma_sqrt[,Ajhat])%*%((1/2)*((as.matrix(B_true[,i]))%*%(B_true[jAjhat,])+(as.matrix(B_true[,jAjhat]))%*%(B_true[i,])))%*%Sigma_sqrt[Ajhat,]
    sigmahat_ij <- sqrt((((B[i,])%*%t(as.matrix(X_plus[,Ajhat])))^2%*%((as.matrix(X_plus[,Ajhat]))%*%(as.matrix(B[,jAjhat])))^2)/n-(1/2)*(theta_hat_ij^2+theta_hat_ji^2))
    sigmatrue_ij <- sqrt(((2*(sqrt(3))^4)/9)*(sum(diag(A%*%A)))-((2*(sqrt(3))^4)/15)*(sum((diag(A))^2)))
    z_ij <- sqrt(n)*(theta_hat_ij-omega2[,j][Ajhat][i])/sigmahat_ij
    z_unbias <- c(z_unbias,z_ij)
    l_ij <- 2*qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
    l_unbias <- c(l_unbias,l_ij)
    ltrue_ij <- 2*qnorm(1-alpha/2)*sigmatrue_ij/sqrt(n)
    ltrue_unbias <- c(ltrue_unbias,ltrue_ij)
    CI_ij_low <- theta_hat_ij-qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
    CI_ij_up <- theta_hat_ij+qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
    if((omega2[,j][Ajhat][i])<=CI_ij_up & (omega2[,j][Ajhat][i])>=CI_ij_low){
      cov_ij <- 1
    }
    else{cov_ij <- 0}
    cov_unbias <- c(cov_unbias,cov_ij)
  }
}

save(Ahat,cov_desparsified,z_desparsified,l_desparsified,ltrue_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("L0_desparsified_unbias_subGaussian_band_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov_desparsified,z_desparsified,l_desparsified,ltrue_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("L0_desparsified_unbias_subGaussian_random_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov_desparsified,z_desparsified,l_desparsified,ltrue_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("L0_desparsified_unbias_subGaussian_hub_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov_desparsified,z_desparsified,l_desparsified,ltrue_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("L0_desparsified_unbias_subGaussian_cluster_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))


