#Run L0BnB algorithm
####################################
#X-----------Design matrix(centralized)
#y-----------Response variable
#kmax--------Max number of nonzero elements in coefficient
#n-----------Sample size
#p-----------Dimension
#m_multi-----A scalar >= 1 used in estimating the big-M. The big-M is defined as the L_infinity norm of the warm start times m_mulitplier. Defaults to 1.2. Larger values can increase the run time.
#time_lim----Maximum run time in seconds per solution.
#tol---------Tolerance
#lambda0r----A multiplier to update the lambda0
#max_iter----To terminate the algorithm early if the number of iterations reaches max_iter
####################################
#A function returns the optimal estimation of coefficient using HBIC criterion
abnb<-function(X,y,kmax,n,p,m_multi,time_lim,tol,lambda0r,max_iter){
  X<-as.matrix(X)
  y<-as.vector(y)
  sols <- fit_path(X, y, lambda_2 = 0, max_nonzeros = kmax, intercept=FALSE, m_multiplier = m_multi, time_limit=time_lim, gap_tol = tol, lambda_0_ratio=lambda0r, max_iteration=max_iter)
  k_grid <- vector()
  k_length <- length(sols)
  for(i in 1:k_length){
    k_grid[i] <- sum(sols[[i]]$B!=0)
  }
  
  HBic={}
  k=0
  Newbeta<-matrix(0,p-1,1)
  #Using HBIC criterion tunning parameters
  while(k<=k_length){
    #k=0 corresponds to null model
    if(k==0){
      ebeta <- as.matrix(rep(0,p-1),ncol=1)
      Newbeta<-cbind(Newbeta,ebeta)
      residual<-(base::norm(X%*%ebeta-y,"2"))^2/n
      hbick<-log(residual)+log(log(n))*log(p-1)*k/n  
      HBic<-rbind(HBic,hbick) 
    }
    else{
      ebeta<-as.vector(sols[[k]]$B)
      Newbeta<-cbind(Newbeta,ebeta)
      residual<-(base::norm(X%*%ebeta-y,"2"))^2/n
      hbick<-log(residual)+log(log(n))*log(p-1)*k_grid[k]/n  
      HBic<-rbind(HBic,hbick)
    }
    k=k+1
  }
  Newbeta<-Newbeta[,-1]
  return(Newbeta[,which.min(HBic)])
}