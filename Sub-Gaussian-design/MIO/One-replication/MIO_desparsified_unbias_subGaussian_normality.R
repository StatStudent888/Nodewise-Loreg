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
Beta <- {}
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
  Beta<-cbind(Beta,best_ebeta)
}
Ahat <- which(Beta!=0)
sigma<-sigma[-1]
sigma_c<-sigma_c[,-1]
#Omega_hat_US
MIO_F_sigma<-matrix_generator3(sigma,sigma_c) #No symmetrization
#T_hat
TT <- matrixMultiply(t(MIO_F_sigma),X_plus2)
MIO_T<-MIO_F_sigma+t(MIO_F_sigma)-matrixMultiply(TT,MIO_F_sigma)

cov_desparsified <- {} #Coverage rate
z_desparsified <- {} #Z-score
l_desparsified <- {} #Avglength
ltrue_desparsified <- {} #Truelength
for(j in 1:p){
  for(i in 1:p){
    T_hat_ij <- MIO_T[i,j]
    A <- Sigma_sqrt%*%((1/2)*(as.matrix(omega2[,i])%*%omega2[j,]+as.matrix(omega2[,j])%*%omega2[i,]))%*%Sigma_sqrt
    sigmatrue_ij <- sqrt(((2*(sqrt(3))^4)/9)*(sum(diag(A%*%A)))-((2*(sqrt(3))^4)/15)*(sum((diag(A))^2)))
    ltrue_ij <- 2*qnorm(1-alpha/2)*sigmatrue_ij/sqrt(n)
    ltrue_desparsified <- c(ltrue_desparsified,ltrue_ij)
    
    sigmahat_ij_lowertri2_left <- ((t(as.matrix(MIO_F_sigma[,i]))%*%t(data))^2%*%(data%*%as.matrix(MIO_F_sigma[,j]))^2)/n
    sigmahat_ij_lowertri2 <- sigmahat_ij_lowertri2_left-(1/2)*((MIO_F_sigma[i,j])^2+(MIO_F_sigma[j,i])^2)
    if(sigmahat_ij_lowertri2<=0){
      sigmahat_ij_lowertri2_fix <- sigmahat_ij_lowertri2_left-min((MIO_F_sigma[i,j])^2,(MIO_F_sigma[j,i])^2)
      if(sigmahat_ij_lowertri2_fix<=0){
        z_ij <- NA
        z_desparsified <- c(z_desparsified,z_ij)
        l_ij <- NA
        l_desparsified <- c(l_desparsified,l_ij)
        cov_ij <- NA
        cov_desparsified <- c(cov_desparsified,cov_ij)
      }
      else{
        sigmahat_ij <- sqrt(sigmahat_ij_lowertri2_fix)
        z_ij <- sqrt(n)*(T_hat_ij-omega2[i,j])/sigmahat_ij
        z_desparsified <- c(z_desparsified,z_ij)
        l_ij <- 2*qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
        l_desparsified <- c(l_desparsified,l_ij)
        
        CI_ij_low <- T_hat_ij-qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
        CI_ij_up <- T_hat_ij+qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
        if((omega2[i,j])<=CI_ij_up & (omega2[i,j])>=CI_ij_low){
          cov_ij <- 1
        }
        else{cov_ij <- 0}
        cov_desparsified <- c(cov_desparsified,cov_ij)
      }
    }
    else{
      sigmahat_ij <- sqrt(sigmahat_ij_lowertri2)
      z_ij <- sqrt(n)*(T_hat_ij-omega2[i,j])/sigmahat_ij
      z_desparsified <- c(z_desparsified,z_ij)
      l_ij <- 2*qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
      l_desparsified <- c(l_desparsified,l_ij)
      
      CI_ij_low <- T_hat_ij-qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
      CI_ij_up <- T_hat_ij+qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
      if((omega2[i,j])<=CI_ij_up & (omega2[i,j])>=CI_ij_low){
        cov_ij <- 1
      }
      else{cov_ij <- 0}
      cov_desparsified <- c(cov_desparsified,cov_ij)
    }
  }
}

cov_unbias <- {}
z_unbias <- {}
l_unbias <- {}
ltrue_unbias <- {}
j <- 0
while(j<p){
  j<-j+1
  Ajhat <- which(MIO_F_sigma[,j]!=0)
  nj <- length(Ajhat)
  jAjhat <- which(Ajhat==j)
  if(nj==0){next}
  for(i in 1:nj){
    theta_hat_ij <- MIO_F_sigma[,j][Ajhat][i]
    theta_hat_ji <- t(MIO_F_sigma)[,j][Ajhat][i]
    B <- solve(X_plus2[Ajhat,Ajhat])
    B_true <- solve(Sigma2[Ajhat,Ajhat])
    A <- as.matrix(Sigma_sqrt[,Ajhat])%*%((1/2)*((as.matrix(B_true[,i]))%*%(B_true[jAjhat,])+(as.matrix(B_true[,jAjhat]))%*%(B_true[i,])))%*%Sigma_sqrt[Ajhat,]
    sigmatrue_ij <- sqrt(((2*(sqrt(3))^4)/9)*(sum(A^2))-((2*(sqrt(3))^4)/15)*(sum((diag(A))^2)))
    ltrue_ij <- 2*qnorm(1-alpha/2)*sigmatrue_ij/sqrt(n)
    ltrue_unbias <- c(ltrue_unbias,ltrue_ij)
    
    sigmahat_ij_lowertri2_left <- (((B[i,])%*%t(as.matrix(X_plus[,Ajhat])))^2%*%(as.matrix(X_plus[,Ajhat])%*%(as.matrix(B[,jAjhat])))^2)/n
    sigmahat_ij_lowertri2 <- sigmahat_ij_lowertri2_left-(1/2)*(theta_hat_ij^2+theta_hat_ji^2)
    if(sigmahat_ij_lowertri2<=0){
      sigmahat_ij_lowertri2_fix <- sigmahat_ij_lowertri2_left-min(theta_hat_ij^2,theta_hat_ji^2)
      if(sigmahat_ij_lowertri2_fix<=0){
        z_ij <- NA
        z_unbias <- c(z_unbias,z_ij)
        l_ij <- NA
        l_unbias <- c(l_unbias,l_ij)
        cov_ij <- NA
        cov_unbias <- c(cov_unbias,cov_ij)
      }
      else{
        sigmahat_ij <- sqrt(sigmahat_ij_lowertri2_fix)
        z_ij <- sqrt(n)*(theta_hat_ij-omega2[,j][Ajhat][i])/sigmahat_ij
        z_unbias <- c(z_unbias,z_ij)
        l_ij <- 2*qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
        l_unbias <- c(l_unbias,l_ij)
        CI_ij_low <- theta_hat_ij-qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
        CI_ij_up <- theta_hat_ij+qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
        if((omega2[,j][Ajhat][i])<=CI_ij_up & (omega2[,j][Ajhat][i])>=CI_ij_low){
          cov_ij <- 1
        }
        else{cov_ij <- 0}
        cov_unbias <- c(cov_unbias,cov_ij)
      }
    }
    else{
      sigmahat_ij <- sqrt(sigmahat_ij_lowertri2)
      z_ij <- sqrt(n)*(theta_hat_ij-omega2[,j][Ajhat][i])/sigmahat_ij
      z_unbias <- c(z_unbias,z_ij)
      l_ij <- 2*qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
      l_unbias <- c(l_unbias,l_ij)
      CI_ij_low <- theta_hat_ij-qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
      CI_ij_up <- theta_hat_ij+qnorm(1-alpha/2)*sigmahat_ij/sqrt(n)
      if((omega2[,j][Ajhat][i])<=CI_ij_up & (omega2[,j][Ajhat][i])>=CI_ij_low){
        cov_ij <- 1
      }
      else{cov_ij <- 0}
      cov_unbias <- c(cov_unbias,cov_ij)
    }
  }
}

save(Ahat,cov_desparsified,z_desparsified,l_desparsified,ltrue_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("MIO_desparsified_unbias_subGaussian_band_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov_desparsified,z_desparsified,l_desparsified,ltrue_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("MIO_desparsified_unbias_subGaussian_random_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov_desparsified,z_desparsified,l_desparsified,ltrue_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("MIO_desparsified_unbias_subGaussian_hub_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))

#save(Ahat,cov_desparsified,z_desparsified,l_desparsified,ltrue_desparsified,cov_unbias,z_unbias,l_unbias,ltrue_unbias,file=paste("MIO_desparsified_unbias_subGaussian_cluster_n",n,"_p",p,"_seed",seed,"_normality.RData",sep=""))






