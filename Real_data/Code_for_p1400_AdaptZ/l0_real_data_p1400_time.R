################################# 
##Real data
#################################
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(caret)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("use_for_all/asdar_youhua_Eigen.cpp")
source("use_for_all/matrix_generator.R")
source("use_for_all/score.R")
source("use_for_all/real_evaluation.R")
source("use_for_all/useful_function.R")
source("use_for_all/FDR-R-Code/epsest.func.R")
source("use_for_all/FDR-R-Code/lin.itp.R")
source("use_for_all/FDR-R-Code/adpt.cutz.R")
source("use_for_all/FDR-R-Code/adaptZ.func.R")

#Read data 
##############################################################
PCR_log2_genedata<-read.csv("Log2/PCR_Test_Log2_genedata_lfdr1400.csv",sep = ",")
RD_log2_genedata<-read.csv("Log2/RD_Test_Log2_genedata_lfdr1400.csv",sep = ",")
PCR_log2_genedata<-PCR_log2_genedata[,-1]
RD_log2_genedata<-RD_log2_genedata[,-1]

#FDR level
alpha <- 0.05
#Dimension
p <- 1400

args <- commandArgs(TRUE)
seed_id=strtoi(args[1])
SEED<-scan("realdata_seed.txt")
seed = SEED[1] #1-100
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

resultreal <- normalizeMatrix(data)
data_normalize <- resultreal$normalizedX
diag <- resultreal$norms

time<-system.time({
  result <- process(data,data_normalize,diag,floor(n/(log(p)*log(log(n)))),n,p)
  sdar_F_sigma1<-matrix_generator1(as.vector(result$sigma),result$sigma_c) }) 

########Compute LDA score on test set########################################################################################
K<-nrow(test)
k<-0
#LDA score for Omega_hat_S
derta_1_min<-matrix()
derta_2_min<-matrix()
while (k<K) {
  k<-k+1
  x<-test[k,]
  tx<-t(x)
  derta_k_min<-score(x,sdar_F_sigma1,mean_PCR,mean_RD,PIP,PIR)
  derta_1_min<- rbind(derta_1_min,derta_k_min$derta1)
  derta_2_min<- rbind(derta_2_min,derta_k_min$derta2)
}
derta_min<-cbind(derta_1_min,derta_2_min)
derta_min<-derta_min[-1,]
derta_pre_min<-matrix()  
k<-0
while(k<K){
  k<-k+1
  derta_prek_min<-which.max(derta_min[k,])
  derta_pre_min<-rbind(derta_pre_min,derta_prek_min)
}

prek_min<-which(derta_pre_min==2)
derta_pre_min[prek_min,]<-0
derta_pre_min<-derta_pre_min[-1,]

derta_true<-matrix(0,K,1)
derta_true[1:5,]<-1  
#Compute classification performance
eva_min <- real_evaluation(derta_pre_min,derta_true)
sparsity <- length(which(sdar_F_sigma1[upper.tri(sdar_F_sigma1)]!=0))

sdar_time<-time
precision <- eva_min$precision
specificity<-eva_min$specificity
sensitivity<-eva_min$sensitivity
MCC<-eva_min$MCC

save(sdar_time,sparsity,precision,specificity,sensitivity,MCC,file=paste("L0_real_data_p",p,"_seed",seed,".RData",sep=""))
