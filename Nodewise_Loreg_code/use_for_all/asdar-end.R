#Run SDAR algorithm
#There is no need to input the normalized version of design matrix X, the code will automatically normalize X
####################################
#X-----------Design matrix(centralized)
#y-----------Response variable
#kmax--------Max number of nonzero elements in coefficient
#n-----------Sample size
#p-----------Dimension
##################################
#A function returns the optimal estimation of coefficient using HBIC criterion
asdar<-function(X,y,kmax,n,p){
  X<-as.matrix(X)
  y<-as.matrix(y)
  dx<-vector()
  #Normalized
  nX<-X
  for (i in 1:ncol(X)){
    dx[i]<-1/norm(X[,i],"2")
    nX[,i]<-X[,i]*dx[i]
  }
  D<-diag(dx)
  HBic={}
  k=0
  Newbeta<-matrix(0,p-1,1)
  #Using HBIC criterion tunning parameters
  while(k<=kmax){
    #k=0 corresponds to null model
    if(k==0){
      ebeta <- as.matrix(rep(0,p-1),ncol=1)
      Newbeta<-cbind(Newbeta,ebeta)
      residual<-(norm(X%*%ebeta-y,"2"))^2/n
      hbick<-log(residual)+log(log(n))*log(p-1)*k/n  
      HBic<-rbind(HBic,hbick) 
    }
    else{
      a<-sdar(nX,y,k,20,0) 
      ebeta<-D%*%a[[3]]
      Newbeta<-cbind(Newbeta,ebeta)
      residual<-(norm(X%*%ebeta-y,"2"))^2/n
      hbick<-log(residual)+log(log(n))*log(p-1)*k/n  
      HBic<-rbind(HBic,hbick)
    }
    k=k+1
  }
  Newbeta<-Newbeta[,-1]
  return(Newbeta[,which.min(HBic)])
}

