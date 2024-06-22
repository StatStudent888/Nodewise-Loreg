#Run Lasso algorithm
####################################
#X-----------Design matrix(centralized)
#y-----------Response variable
#kmax--------Max number of nonzero elements in coefficient
#n-----------Sample size
#p-----------Dimension
#bestbeta----Optimal coefficient using HBIC criterion
#bestlambda--Optimal parameter using HBIC criterion
####################################
#A function returns the optimal estimation of coefficient using HBIC criterion
ans<-function(X,y,Lambda,n,p){
  X<-as.matrix(X)
  y<-as.matrix(y)
  #nX<-scale(X)
  HBic={}
  betaold=matrix(0,p-1,1)  
  k=0
  Newbeta<-matrix(0,p-1,1)  
  Newlambda <- {}
  kmax<-length(Lambda)
  M<-{}
  #Using HBIC criterion tunning parameters
  while(k<kmax){
    k=k+1
    l=Lambda[k]
    Newlambda <- c(Newlambda,l)
    model<-glmnet(X,y,family="gaussian",alpha=1,lambda =l)
    ebeta = as.matrix(model$beta)
    residual<-(norm(X%*%ebeta-y,"2"))^2/n
    Newbeta<-cbind(Newbeta,ebeta)
    m<-length(which(ebeta!=0))
    hbick<-log(residual) +log(log(n))*log(p-1)*m/n  
    HBic<-rbind(HBic,hbick)         
  }
  Newbeta=Newbeta[,-1]  
  bestbeta=Newbeta[,which.min(HBic)]
  bestlambda=Newlambda[which.min(HBic)]
  return(list(bestbeta=bestbeta,bestlambda=bestlambda))
}

