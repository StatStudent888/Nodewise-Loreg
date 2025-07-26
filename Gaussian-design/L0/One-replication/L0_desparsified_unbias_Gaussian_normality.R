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

SEED<-scan("seed.txt")
seed = SEED[1] #1-100
set.seed(seed)
#Generate Gaussian data
data <- mvrnorm(n,Mean,Sigma2)
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
TT <- matrixMultiply(t(sdar_F_sigma),X_plus2)
#T_hat
sdar_T<-sdar_F_sigma+t(sdar_F_sigma)-matrixMultiply(TT,sdar_F_sigma) #Desparsified estimator

result <- desparsified_confidence_intervals_Gaussian(sdar_T,sdar_F_sigma,omega2,n,p,alpha)
cov_desparsified <- result$cov_desparsified #Coverage rate
z_desparsified <- result$z_desparsified #Z-score
l_desparsified <- result$l_desparsified #Avglength

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
    sigmahat_ij <- sqrt((solve(X_plus2[Ajhat,Ajhat])[i,i])*(solve(X_plus2[Ajhat,Ajhat])[jAjhat,jAjhat])+(solve(X_plus2[Ajhat,Ajhat])[i,jAjhat])^2)
    sigmatrue_ij <- sqrt((solve(Sigma2[Ajhat,Ajhat])[i,i])*(solve(Sigma2[Ajhat,Ajhat])[jAjhat,jAjhat])+(solve(Sigma2[Ajhat,Ajhat])[i,jAjhat])^2)
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

z_desparsified <- round(z_desparsified,3)
l_desparsified <- round(l_desparsified,3)
z_unbias <- round(z_unbias,3)
l_unbias <- round(l_unbias,3)
ltrue_unbias <- round(ltrue_unbias,3)

save(Ahat,cov_desparsified,z_desparsified,l_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("L0_desparsified_unbias_Gaussian_band_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

# save(Ahat,cov_desparsified,z_desparsified,l_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("L0_desparsified_unbias_Gaussian_random_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))
# 
# save(Ahat,cov_desparsified,z_desparsified,l_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("L0_desparsified_unbias_Gaussian_hub_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))
# 
# save(Ahat,cov_desparsified,z_desparsified,l_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("L0_desparsified_unbias_Gaussian_cluster_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

