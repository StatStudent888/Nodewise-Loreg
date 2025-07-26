library(glasso)
library(huge)
library(MASS)
source("use_for_all/omega_generator.R")
source("use_for_all/evaluation.R")

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

#Generate Sub-Gaussian data
Sigma_sqrt <- matrix_sqrt(Sigma)
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
X<-data
X<-as.matrix(X)

Bestwi<-list()
BIC<-{}

#Estimate precision matrix and compute running time
timek<-system.time({
  S<-var(X)  
  for(j in 1:20){
    l=Lambda[j]
    fit<-glasso(S,rho=l)
    fit_wi<-fit$wi
    Bestwi[[j]]<-fit_wi
    C<-as.matrix(fit_wi)  
    c=matrix(0,nrow=p,ncol=p)
    c[row(c)<=col(c)]=C[row(C)<=col(C)]
    bick<-sum(diag(C%*%S))-as.numeric(determinant(C)$modulus)+(log(n)/n)*(sum(c!=0))  
    BIC<- c(BIC,bick) 
  }
  #******************************
  ##%% Output  a solution 
  id<-which.min(BIC)
  Glasso_F_sigma<- Bestwi[[id]]
  Glasso_F_sigma<- as.matrix(Glasso_F_sigma)}) 

#Support recovery performance
eva <- evaluation(omega,Glasso_F_sigma)
precision <- eva$precision
sensitivity <- eva$sensitivity
specificity <- eva$specificity
MCC <- eva$MCC

#Matrix norm loss
Glasso_norm1<-norm(Glasso_F_sigma-omega,"1")
Glasso_norm2<-norm(Glasso_F_sigma-omega,"2")
Glasso_normF<-norm(Glasso_F_sigma-omega,"F")
Glasso_normmax<-norm(Glasso_F_sigma-omega,"M")

#Running time
Glasso_time<-timek

save(precision,sensitivity,specificity,MCC,Glasso_norm1,Glasso_norm2,Glasso_normF,Glasso_normmax,Glasso_time,file=paste("Glasso_subGaussian_band_n",n,"_p",p,"_seed",seed,".RData",sep=""))

# save(precision,sensitivity,specificity,MCC,Glasso_norm1,Glasso_norm2,Glasso_normF,Glasso_normmax,Glasso_time,file=paste("Glasso_subGaussian_random_n",n,"_p",p,"_seed",seed,".RData",sep=""))
# 
# save(precision,sensitivity,specificity,MCC,Glasso_norm1,Glasso_norm2,Glasso_normF,Glasso_normmax,Glasso_time,file=paste("Glasso_subGaussian_hub_n",n,"_p",p,"_seed",seed,".RData",sep=""))
# 
# save(precision,sensitivity,specificity,MCC,Glasso_norm1,Glasso_norm2,Glasso_normF,Glasso_normmax,Glasso_time,file=paste("Glasso_subGaussian_cluster_n",n,"_p",p,"_seed",seed,".RData",sep=""))






