#Run MIO algorithm
####################################
#X-----------Design matrix(centralized)
#y-----------Response variable
#kmax--------Max number of nonzero elements in coefficient
#n-----------Sample size
#p-----------Dimension
#type_penal--Type of penalty: l1 or l2, which is useless when we set llambda=0, the corresponding coefficient for lq penalty is zero
#time_lim----Maximum run time in seconds per solution.
####################################
#A function returns the optimal estimation of coefficient using HBIC criterion
amio<-function(X,y,kmax,n,p,type_penal,time_lim){
  X<-as.matrix(X)
  y<-as.vector(y)
  dx<-vector()
  #Normalized
  nX<-X
  for (i in 1:ncol(X)){
    dx[i]<-1/base::norm(X[,i],"2")
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
      residual<-(base::norm(X%*%ebeta-y,"2"))^2/n
      hbick<-log(residual)+log(log(n))*log(p-1)*k/n  
      HBic<-rbind(HBic,hbick) 
    }
    else{
      a<-MIO_L0_LQ(type_penal, nX, y, r_to_py(as.integer(k)), 0, time_limit = time_lim) #llambda=0, the coefficient for the lq penalty
      ebeta<-as.vector(a[[1]]%*%D)
      Newbeta<-cbind(Newbeta,ebeta)
      residual<-(base::norm(X%*%ebeta-y,"2"))^2/n
      hbick<-log(residual)+log(log(n))*log(p-1)*k/n  
      HBic<-rbind(HBic,hbick)
    }
    k=k+1
  }
  Newbeta<-Newbeta[,-1]
  return(Newbeta[,which.min(HBic)])
}