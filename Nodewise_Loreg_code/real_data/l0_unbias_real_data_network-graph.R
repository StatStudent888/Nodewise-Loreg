################################# 
##Real data network graph
#################################
library(Rcpp)
library(RcppArmadillo)
library(caret)
library(hgu133plus2.db)
library(ggraph)
library(igraph)
library(extras)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("end-sdar/sdar_RCPP-end.cpp")
source("use_for_all/asdar-end.R")
source("use_for_all/matrix_generator.R")
source("FDR-R-Code/epsest.func.R")
source("FDR-R-Code/lin.itp.R")
source("FDR-R-Code/adpt.cutz.R")
source("FDR-R-Code/adaptZ.func.R")
#Read data 
##############################################################
PCR_log2_genedata<-read.csv("Log2/PCR_Test_Log2_genedata_300.csv",sep = ",")
RD_log2_genedata<-read.csv("Log2/RD_Test_Log2_genedata_300.csv",sep = ",")
probe_id <- PCR_log2_genedata[,1]
PCR_log2_genedata<-PCR_log2_genedata[,-1]
RD_log2_genedata<-RD_log2_genedata[,-1]
NP<-ncol(PCR_log2_genedata)
NR<-ncol(RD_log2_genedata)

#FDR level
alpha <- 0.05
#Dimension
p <- 300

mean_PCR<-matrix(0,p,1)
mean_RD<-matrix(0,p,1)
M_PCR_log2_genedata<-PCR_log2_genedata
M_RD_log2_genedata<-RD_log2_genedata

for(ei in 1:p){
  MP<-sum(PCR_log2_genedata[ei,])/NP
  MR<-sum(RD_log2_genedata[ei,])/NR
  mean_PCR[ei,]<-MP
  mean_RD[ei,]<-MR
  M_PCR_log2_genedata[ei,]<-PCR_log2_genedata[ei,]-MP
  M_RD_log2_genedata[ei,]<-RD_log2_genedata[ei,]-MR
}

data <- t(cbind(M_PCR_log2_genedata,M_RD_log2_genedata))
data<-as.matrix(data)
n <- nrow(data)

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
  best_ebeta<-asdar(X,y,14,n,p) #14=floor(n/(logploglogn)),this time we use whole dataset
  sigmak<-n/(norm(y-X%*%best_ebeta,"2")^2)
  sigma_ck<-(-1)*sigmak*best_ebeta
  sigma<- cbind(sigma,sigmak)
  sigma_c<- cbind(sigma_c,sigma_ck)
}
sigma<-sigma[-1]
sigma_c<-sigma_c[,-1]
#Omega_hat_S
sdar_F_sigma1<-matrix_generator1(sigma,sigma_c)
#Omega_hat_US
sdar_F_sigma3<-matrix_generator3(sigma,sigma_c)
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
sdar_unbias_F_sigma_star <- matrix_generator3(sigma,sigma_c_star1) 

colnames(sdar_unbias_F_sigma_star) <- probe_id
rownames(sdar_unbias_F_sigma_star) <- probe_id

#Compute edges of T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S))
sparsity_unbias_s_min <- length(which(sdar_unbias_F_sigma_star[upper.tri(sdar_unbias_F_sigma_star)]!=0))
sparsity_unbias_s_min

#Estimated absolute partial correlation matrix
sdar_unbias_partial_correlation_star <- matrix(0,nrow=p,ncol=p)
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){
      sdar_unbias_partial_correlation_star[i,j] <- abs(-sdar_unbias_F_sigma_star[i,j]/(sqrt(sdar_unbias_F_sigma_star[i,i]*sdar_unbias_F_sigma_star[j,j])))
    }
    else{sdar_unbias_partial_correlation_star[i,j] <- 0}
  }
}
colnames(sdar_unbias_partial_correlation_star) <- probe_id
rownames(sdar_unbias_partial_correlation_star) <- probe_id

#Network graph of estimated absolute partial correlation matrix
library(igraph)
g1 <- graph_from_adjacency_matrix(sdar_unbias_partial_correlation_star,mode = "undirected",weighted=TRUE)
V(g1)$size=igraph::degree(g1)%>%+1%>%pow(0.35)*2
E(g1)$width=E(g1)$weight%>%+1%>%pow(2)
set.seed(1)
l=layout_nicely(g1)
plot(g1, vertex.label.color="cornflowerblue", 
     vertex.label.cex=0.6, vertex.label.dist=0.5,vertex.frame.color="azure4",layout=l)

#Estimated absolute partial correlation matrix thresholded by 0.1
sdar_partial_correlation_threshold1 <- sdar_unbias_partial_correlation_star
sdar_partial_correlation_threshold1[which(sdar_unbias_partial_correlation_star<=0.1)] <- 0
#Compute edges of estimated absolute partial correlation matrix thresholded by 0.1
sparsity_1 <- length(which(sdar_partial_correlation_threshold1[upper.tri(sdar_partial_correlation_threshold1)]!=0))
sparsity_1

#Network graph of estimated absolute partial correlation matrix thresholded by 0.1
g1_t1 <- graph_from_adjacency_matrix(sdar_partial_correlation_threshold1,mode = "undirected",weighted=TRUE)
V(g1_t1)$size=igraph::degree(g1_t1)%>%+1%>%pow(0.35)*2
E(g1_t1)$width=E(g1_t1)$weight%>%+1%>%pow(2)
plot(g1_t1, vertex.label.color="cornflowerblue", 
     vertex.label.cex=0.6, vertex.label.dist=0.5,vertex.frame.color="azure4",layout=l)


