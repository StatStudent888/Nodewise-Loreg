library(glmnet)
library(huge)
library(MASS)
source("use_for_all/ans-end.R")
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

#Candidate values of tunning parameter
num_lambda=20;lambda_max=2;b=50;lambda_min_ratio=0.01
lambdamax = lambda_max
nlambda=num_lambda
lambda.min.ratio=lambda_min_ratio
lambda.min = lambda.min.ratio*lambdamax
Lambda = exp(seq(log(lambdamax), log(lambda.min), length = nlambda))

SEED<-scan("seed.txt")
seed = SEED[1]
set.seed(seed)
#Generate Gaussian data
data_Gaussian <- mvrnorm(n,Mean,Sigma2)

#Generate sub-Gaussian data
Sigma_sqrt <- eigen(Sigma2)$vectors%*%sqrt(diag(eigen(Sigma2)$values))%*%solve(eigen(Sigma2)$vectors)
U <- matrix(0,n,p)
set.seed(seed)
for(h in 1:n){
  for(q in 1:p){
    U[h,q] <- runif(1,-sqrt(3),sqrt(3))
  }
}
data_subGaussian <- U%*%Sigma_sqrt

j<-0
sigma_Gaussian<-matrix()
sigma_c_Gaussian<-matrix(1,p-1,1)

#Time for Gaussian setting
timek_Gaussian<-system.time({
  while(j<p){
    j<-j+1
    X<-data_Gaussian[,-j]
    y<-data_Gaussian[,j]
    X<-as.matrix(X) 
    y<-as.vector(y)
    best_ebeta_lambda<-ans(X,y,Lambda,n,p)
    best_ebeta <- best_ebeta_lambda$bestbeta
    best_lambda <- best_ebeta_lambda$bestlambda
    sigmak<-1/(((norm(y-X%*%best_ebeta,"2")^2)/n)+best_lambda/2*sum(abs(best_ebeta))) 
    sigma_ck<-(-1)*sigmak*best_ebeta  
    sigma_Gaussian<- cbind(sigma_Gaussian,sigmak)
    sigma_c_Gaussian<- cbind(sigma_c_Gaussian,sigma_ck)
  }
  sigma_Gaussian<-sigma_Gaussian[-1]
  sigma_c_Gaussian<-sigma_c_Gaussian[,-1]
  NS_F_sigma1_Gaussian<-matrix_generator1(sigma_Gaussian,sigma_c_Gaussian) }) 

j<-0
sigma_subGaussian<-matrix()
sigma_c_subGaussian<-matrix(1,p-1,1)

#Time for sub-Gaussian setting
timek_subGaussian<-system.time({
  while(j<p){
    j<-j+1
    X<-data_subGaussian[,-j]
    y<-data_subGaussian[,j]
    X<-as.matrix(X) 
    y<-as.vector(y)
    best_ebeta_lambda<-ans(X,y,Lambda,n,p)
    best_ebeta <- best_ebeta_lambda$bestbeta
    best_lambda <- best_ebeta_lambda$bestlambda
    sigmak<-1/(((norm(y-X%*%best_ebeta,"2")^2)/n)+best_lambda/2*sum(abs(best_ebeta))) 
    sigma_ck<-(-1)*sigmak*best_ebeta  
    sigma_subGaussian<- cbind(sigma_subGaussian,sigmak)
    sigma_c_subGaussian<- cbind(sigma_c_subGaussian,sigma_ck)
  }
  sigma_subGaussian<-sigma_subGaussian[-1]
  sigma_c_subGaussian<-sigma_c_subGaussian[,-1]
  NS_F_sigma1_subGaussian<-matrix_generator1(sigma_subGaussian,sigma_c_subGaussian) })


save(data_Gaussian,data_subGaussian,NS_F_sigma1_Gaussian,NS_F_sigma1_subGaussian,timek_Gaussian,timek_subGaussian,file=paste("L1_band_n",n,"_p",p,"_seed",seed,"_time.RData",sep=""))

#save(data_Gaussian,data_subGaussian,NS_F_sigma1_Gaussian,NS_F_sigma1_subGaussian,timek_Gaussian,timek_subGaussian,file=paste("L1_random_n",n,"_p",p,"_seed",seed,"_time.RData",sep=""))

#save(data_Gaussian,data_subGaussian,NS_F_sigma1_Gaussian,NS_F_sigma1_subGaussian,timek_Gaussian,timek_subGaussian,file=paste("L1_hub_n",n,"_p",p,"_seed",seed,"_time.RData",sep=""))

#save(data_Gaussian,data_subGaussian,NS_F_sigma1_Gaussian,NS_F_sigma1_subGaussian,timek_Gaussian,timek_subGaussian,file=paste("L1_cluster_n",n,"_p",p,"_seed",seed,"_time.RData",sep=""))