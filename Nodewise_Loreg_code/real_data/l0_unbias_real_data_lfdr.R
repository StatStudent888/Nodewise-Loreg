################################# 
##Real data
#################################
library(Rcpp)
library(RcppArmadillo)
library(caret)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("end-sdar/sdar_RCPP-end.cpp")
source("use_for_all/asdar-end.R")
source("use_for_all/matrix_generator.R")
source("use_for_all/score.R")
source("use_for_all/real_evaluation.R")
source("FDR-R-Code/epsest.func.R")
source("FDR-R-Code/lin.itp.R")
source("FDR-R-Code/adpt.cutz.R")
source("FDR-R-Code/adaptZ.func.R")

#Read data 
##############################################################
PCR_log2_genedata<-read.csv("Log2/PCR_Test_Log2_genedata_300.csv",sep = ",")
RD_log2_genedata<-read.csv("Log2/RD_Test_Log2_genedata_300.csv",sep = ",")
PCR_log2_genedata<-PCR_log2_genedata[,-1]
RD_log2_genedata<-RD_log2_genedata[,-1]

precision_min <- {}
specificity_min<-{}   
sensitivity_min<-{}   
MCC_min<-{}
precision_s_min <- {}
specificity_s_min <- {}
sensitivity_s_min <- {}
MCC_s_min <- {}
sparsity_min <- {}
sparsity_s_min <- {}
#FDR level
alpha <- 0.05
#Dimension
p <- 300

set.seed(-1)
SEED<-runif(100,-100000,100000)
for(t in 1:100){
  seed<-SEED[t]
  set.seed(seed)
  set1<-sample(ncol(PCR_log2_genedata), 5, replace = FALSE, prob =NULL)
  set2<-sample(ncol(RD_log2_genedata), 16, replace = FALSE, prob =NULL)
  #Test set
  test_PCR<-PCR_log2_genedata[,set1]
  test_RD<-RD_log2_genedata[,set2]
  #Train set
  train_PCR<-PCR_log2_genedata[,-set1]
  train_RD<-RD_log2_genedata[,-set2]
  
  #Centralization
  ##########################################
  #Estimate in-class mean based on train set
  #########################################
  NP<-ncol(train_PCR)
  NR<-ncol(train_RD)
  PIP<-NP/(NP+NR)
  PIR<-NR/(NP+NR)
  mean_PCR<-matrix(0,p,1)
  mean_RD<-matrix(0,p,1)
  M_train_PCR<-train_PCR
  M_train_RD<-train_RD
  for(ei in 1:p){
    MP<-sum(train_PCR[ei,])/NP
    MR<-sum(train_RD[ei,])/NR
    mean_PCR[ei,]<-MP
    mean_RD[ei,]<-MR
    #Centralize train set
    M_train_PCR[ei,]<-train_PCR[ei,]-MP
    M_train_RD[ei,]<-train_RD[ei,]-MR
  }
  
  
  ##########################################
  #Combine train set and test set respectively
  #########################################
  train<-t(cbind(M_train_PCR,M_train_RD))
  test<-t(cbind(test_PCR,test_RD))
  
  ##########################################
  #Estimate precision matrix
  #########################################    
  n<-NP+NR
    
  data<-train
  data<-as.matrix(data)
    
  j<-0
  sigma<-matrix()
  sigma_c<-matrix(1,p-1,1)
  while(j<p){
    j<-j+1
    X<-data[,-j]
    y<-data[,j]
    X<-as.matrix(X)
    y<-as.vector(y)
    best_ebeta<-asdar(X,y,12,n,p) #12=floor(n/(logploglogn))
    sigmak<-n/(norm(y-X%*%best_ebeta,"2")^2)
    sigma_ck<-(-1)*sigmak*best_ebeta
    sigma<- cbind(sigma,sigmak)
    sigma_c<- cbind(sigma_c,sigma_ck)
  }
  sigma<-sigma[-1]
  sigma_c<-sigma_c[,-1]
  #Omega_hat_S
  sdar_F_sigma1<-matrix_generator1(sigma,sigma_c)  #Minimum symmetrization
  #Omega_hat_US
  sdar_F_sigma3<-matrix_generator3(sigma,sigma_c) #No symmetrization
  sdar_F_lowertri1 <- sdar_F_sigma1
  sdar_F_lowertri1[upper.tri(sdar_F_lowertri1)] <- 0
  sdar_F_lowertri1_nodiag <- matrix(0,nrow=p-1,ncol=p)
  sdar_F_lowertri1_nodiag[!upper.tri(sdar_F_lowertri1_nodiag)] <- sdar_F_lowertri1[lower.tri(sdar_F_lowertri1)]
  
  #T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S))
  j <- 0
  zvalue <- {}
  while(j<p){
    j<-j+1
    Ajhat <- which(sdar_F_lowertri1[-j,j]!=0) 
    nj <- length(Ajhat)
    Ajhat_plus <- which(sdar_F_sigma3[,j]!=0)
    jAjhat <- which(Ajhat_plus==j)
    if(nj==0){next}
    for(i in 1:nj){
      iAjhat <- which(Ajhat_plus==(Ajhat[i]+1))
      theta_hat_ij <- sdar_F_sigma3[-j,j][Ajhat][i] 
      sigmahat_ij_lowertri <- sqrt((solve((t(data)%*%data/n)[Ajhat_plus,Ajhat_plus])[iAjhat,iAjhat])*(solve((t(data)%*%data/n)[Ajhat_plus,Ajhat_plus])[jAjhat,jAjhat])+(solve((t(data)%*%data/n)[Ajhat_plus,Ajhat_plus])[iAjhat,jAjhat])^2)
      zvalue_ij <- sqrt(n)*theta_hat_ij/sigmahat_ij_lowertri
      zvalue <- c(zvalue,zvalue_ij)
    }
  }
  Ahat <- which(sdar_F_lowertri1_nodiag!=0)
  critical <- Ahat[adaptZ.func(zvalue, alpha)$ac]
  sigma_c_star1 <- sdar_F_lowertri1_nodiag
  sigma_c_star1[critical] <- 0
  sigma_c_star1[upper.tri(sigma_c_star1)] <- t(sigma_c_star1)[!lower.tri(t(sigma_c_star1))]
  #T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S))
  sdar_T_star1 <- matrix_generator3(sigma,sigma_c_star1) 
  
  ########Compute LDA score on test set########################################################################################
  K<-nrow(test)
  k<-0
  #LDA score for Omega_hat_S
  derta_1_min<-matrix()
  derta_2_min<-matrix()
  #T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S))
  derta_1_s_min<-matrix()
  derta_2_s_min<-matrix()
  while (k<K) {
    k<-k+1
    x<-test[k,]
    tx<-t(x)
    derta_k_min<-score(x,sdar_F_sigma1,mean_PCR,mean_RD,PIP,PIR)
    derta_sk_min<-score(x,sdar_T_star1,mean_PCR,mean_RD,PIP,PIR)
    derta_1_min<- rbind(derta_1_min,derta_k_min$derta1)
    derta_2_min<- rbind(derta_2_min,derta_k_min$derta2)
    derta_1_s_min<- rbind(derta_1_s_min,derta_sk_min$derta1)
    derta_2_s_min<- rbind(derta_2_s_min,derta_sk_min$derta2)
  }
  derta_min<-cbind(derta_1_min,derta_2_min)
  derta_s_min<-cbind(derta_1_s_min,derta_2_s_min)
  derta_min<-derta_min[-1,]
  derta_s_min<-derta_s_min[-1,]
  derta_pre_min<-matrix()  
  derta_pre_s_min<-matrix()
  k<-0
  while(k<K){
    k<-k+1
    derta_prek_min<-which.max(derta_min[k,])
    derta_pre_sk_min<-which.max(derta_s_min[k,])
    derta_pre_min<-rbind(derta_pre_min,derta_prek_min)
    derta_pre_s_min<-rbind(derta_pre_s_min,derta_pre_sk_min)
  }
  
  prek_min<-which(derta_pre_min==2)
  pre_sk_min<-which(derta_pre_s_min==2)
  derta_pre_min[prek_min,]<-0
  derta_pre_s_min[pre_sk_min,]<-0
  derta_pre_min<-derta_pre_min[-1,]
  derta_pre_s_min<-derta_pre_s_min[-1,]
  
  derta_true<-matrix(0,K,1)
  derta_true[1:5,]<-1  
  #Compute classification performance
  eva_min <- real_evaluation(derta_pre_min,derta_true)
  eva_s_min <- real_evaluation(derta_pre_s_min,derta_true)
  sparsityk_min <- length(which(sdar_F_sigma1[upper.tri(sdar_F_sigma1)]!=0))
  sparsity_sk_min <- length(which(sdar_T_star1[upper.tri(sdar_T_star1)]!=0))
  
  precision_min <- c(precision_min,eva_min$precision)
  specificity_min<-c(specificity_min,eva_min$specificity)
  sensitivity_min<-c(sensitivity_min,eva_min$sensitivity)
  MCC_min<-c(MCC_min,eva_min$MCC)
  precision_s_min <- c(precision_s_min,eva_s_min$precision)
  specificity_s_min<-c(specificity_s_min,eva_s_min$specificity)
  sensitivity_s_min<-c(sensitivity_s_min,eva_s_min$sensitivity)
  MCC_s_min<-c(MCC_s_min,eva_s_min$MCC)
  sparsity_min <- c(sparsity_min,sparsityk_min)
  sparsity_s_min <- c(sparsity_s_min,sparsity_sk_min)
}
#####Compute mean under 100 replications
M_precision_min <- mean(precision_min,na.rm=T)
M_specificity_min<-mean(specificity_min)
M_sensitivity_min<-mean(sensitivity_min)
M_MCC_min<-mean(MCC_min,na.rm=T)
M_precision_s_min <- mean(precision_s_min)
M_specificity_s_min<-mean(specificity_s_min)
M_sensitivity_s_min<-mean(sensitivity_s_min)
M_MCC_s_min<-mean(MCC_s_min)
M_sparsity_min <- mean(sparsity_min)
M_sparsity_s_min <- mean(sparsity_s_min)

#####Compute sd under 100 replications
sd_precision_min <- sd(precision_min,na.rm=T)
sd_specificity_min<-sd(specificity_min)
sd_sensitivity_min<-sd(sensitivity_min)
sd_MCC_min<-sd(MCC_min,na.rm=T)
sd_precision_s_min <- sd(precision_s_min)
sd_specificity_s_min<-sd(specificity_s_min)
sd_sensitivity_s_min<-sd(sensitivity_s_min)
sd_MCC_s_min<-sd(MCC_s_min)
sd_sparsity_min <- sd(sparsity_min)
sd_sparsity_s_min <- sd(sparsity_s_min)

#Results
M_precision_min 
M_specificity_min
M_sensitivity_min
M_MCC_min
M_precision_s_min 
M_specificity_s_min
M_sensitivity_s_min
M_MCC_s_min
M_sparsity_min 
M_sparsity_s_min 

sd_precision_min 
sd_specificity_min
sd_sensitivity_min
sd_MCC_min
sd_precision_s_min 
sd_specificity_s_min
sd_sensitivity_s_min
sd_MCC_s_min
sd_sparsity_min 
sd_sparsity_s_min 

save.image("L0_unbias_real_data_lfdr_p300.RData")
