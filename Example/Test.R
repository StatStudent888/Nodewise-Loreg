library(Rcpp)
library(RcppArmadillo)
library(huge)
library(MASS)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("use_for_all/asdar.cpp")
source("use_for_all/ans-end.R")
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
source("use_for_all/useful_function.R")

source("use_for_all/FDR-R-Code/epsest.func.R")
source("use_for_all/FDR-R-Code/lin.itp.R")
source("use_for_all/FDR-R-Code/adpt.cutz.R")
source("use_for_all/FDR-R-Code/adaptZ.func.R")

source("Nodewise_Loreg.R")

#Sample size
n=200 
#Dimension
p=50 
#FDR level for multiple test
alpha <- 0.05

#Generate Band graph
set.seed(-1)
omega2 <- omega_generator1(p,1,0.5,0.3)
Sigma2 <- solve(omega2)

Mean <- rep(0,p)

#Generate Gaussian data
set.seed(1)
data <- mvrnorm(n,Mean,Sigma2)

#Candidate values for tunning parameter for Nodewise_Loreg
Tmax <- 20

#Nodewise_Loreg estimator
Nodewise_Loreg_result <- Nodewise_Loregfunc(data,Tmax,alpha)
Nodewise_Loreg_Omega_hat_S <- Nodewise_Loreg_result$Omega_hat_S
Nodewise_Loreg_Omega_hat_star1 <- Nodewise_Loreg_result$Omega_hat_star1
Nodewise_Loreg_Omega_hat_star2 <- Nodewise_Loreg_result$Omega_hat_star2
Nodewise_Loreg_Omega_hat_star3 <- Nodewise_Loreg_result$Omega_hat_star3
Nodewise_Loreg_Omega_hat_star4 <- Nodewise_Loreg_result$Omega_hat_star4
Nodewise_Loreg_T_hat <- Nodewise_Loreg_result$T_hat


