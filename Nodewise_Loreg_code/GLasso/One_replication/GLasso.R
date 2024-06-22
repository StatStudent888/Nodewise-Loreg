library(glasso)
library(huge)
library(MASS)
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")

#Candidate values for tunning parameter
num_lambda=20; lambda_max=2;b=50;lambda_min_ratio=0.01
lambdamax = lambda_max
nlambda=num_lambda
lambda.min.ratio=lambda_min_ratio
lambda.min = lambda.min.ratio*lambdamax
Lambda = exp(seq(log(lambdamax), log(lambda.min), length = nlambda))

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

X_Gaussian<-data_Gaussian
X_Gaussian<-as.matrix(X_Gaussian)
Bestwi_Gaussian<-list()
BIC_Gaussian<-{}

#Time for Gaussian setting
timek_Gaussian<-system.time({
  S_Gaussian<-var(X_Gaussian)  
  for(j in 1:20){  
    l=Lambda[j]
    fit<-glasso(S_Gaussian,rho=l)
    fit_wi<-fit$wi
    Bestwi_Gaussian[[j]]<-fit_wi
    C<-as.matrix(fit_wi)  
    c=matrix(0,nrow=p,ncol=p)
    c[row(c)<=col(c)]=C[row(C)<=col(C)]
    bick<-sum(diag(C%*%S_Gaussian))-log(det(C))+(log(n)/n)*(sum(c!=0))  
    BIC_Gaussian<- c(BIC_Gaussian,bick) 
  }
  #******************************
  ##%% Output  a solution 
  id_Gaussian<-which.min(BIC_Gaussian)
  Glasso_F_sigma_Gaussian<- Bestwi_Gaussian[[id_Gaussian]]
  Glasso_F_sigma_Gaussian<- as.matrix(Glasso_F_sigma_Gaussian)}) 

j<-0
sigma_subGaussian<-matrix()
sigma_c_subGaussian<-matrix(1,p-1,1)

X_subGaussian<-data_subGaussian
X_subGaussian<-as.matrix(X_subGaussian)
Bestwi_subGaussian<-list()
BIC_subGaussian<-{}

#Time for sub-Gaussian setting
timek_subGaussian<-system.time({
  S_subGaussian<-var(X_subGaussian)  
  for(j in 1:20){  
    l=Lambda[j]
    fit<-glasso(S_subGaussian,rho=l)
    fit_wi<-fit$wi
    Bestwi_subGaussian[[j]]<-fit_wi
    C<-as.matrix(fit_wi)  
    c=matrix(0,nrow=p,ncol=p)
    c[row(c)<=col(c)]=C[row(C)<=col(C)]
    bick<-sum(diag(C%*%S_subGaussian))-log(det(C))+(log(n)/n)*(sum(c!=0))  
    BIC_subGaussian<- c(BIC_subGaussian,bick) 
  }
  #******************************
  ##%% Output  a solution 
  id_subGaussian<-which.min(BIC_subGaussian)
  Glasso_F_sigma_subGaussian<- Bestwi_subGaussian[[id_subGaussian]]
  Glasso_F_sigma_subGaussian<- as.matrix(Glasso_F_sigma_subGaussian)}) 


save(data_Gaussian,data_subGaussian,Glasso_F_sigma_Gaussian,Glasso_F_sigma_subGaussian,timek_Gaussian,timek_subGaussian,file=paste("Glasso_band_n",n,"_p",p,"_seed",seed,".RData",sep=""))

#save(data_Gaussian,data_subGaussian,Glasso_F_sigma_Gaussian,Glasso_F_sigma_subGaussian,timek_Gaussian,timek_subGaussian,file=paste("Glasso_random_n",n,"_p",p,"_seed",seed,".RData",sep=""))

#save(data_Gaussian,data_subGaussian,Glasso_F_sigma_Gaussian,Glasso_F_sigma_subGaussian,timek_Gaussian,timek_subGaussian,file=paste("Glasso_hub_n",n,"_p",p,"_seed",seed,".RData",sep=""))

#save(data_Gaussian,data_subGaussian,Glasso_F_sigma_Gaussian,Glasso_F_sigma_subGaussian,timek_Gaussian,timek_subGaussian,file=paste("Glasso_cluster_n",n,"_p",p,"_seed",seed,".RData",sep=""))