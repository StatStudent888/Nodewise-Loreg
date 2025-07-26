#Input
#1. Data matrix: data(sample size n, dimension p)
#2. Candidate value for tuning parameter: Tmax
#3. fdr level for multiple test: alpha

#Output
#Six different Nodewise_Loreg estimators (see our paper: Nodewise Loreg: Nodewise L0-penalized Regression for High-dimensional Sparse Precision Matrix Estimation)
Nodewise_Loregfunc <- function(data,Tmax,alpha){
  n <- nrow(data)
  p <- ncol(data)
  X_plus <- as.matrix(data)
  X_plus2 <- matrixMultiply(t(X_plus),X_plus)/n
  
  result <- normalizeMatrix(data)
  data_normalize <- result$normalizedX
  diag <- result$norms
  
  j<-0
  sigma<-matrix()
  sigma_c<-matrix(1,p-1,1)
  while(j<p){
    j<-j+1
    X<-data[,-j]
    y<-data[,j]
    nX <- data_normalize[,-j]
    dx <- diag[-j]
    y<-as.vector(y)
    best_ebeta<-asdar(X,nX,dx,y,Tmax,n,p)
    sigmak<-n/(norm(y-X%*%best_ebeta,"2")^2)  
    sigma_ck<-(-1)*sigmak*best_ebeta   
    sigma<- cbind(sigma,sigmak)
    sigma_c<- cbind(sigma_c,sigma_ck)
  }
  sigma<-sigma[-1]
  sigma_c<-sigma_c[,-1]
  #Omega_hat_S
  sdar_F_sigma1<-matrix_generator1(sigma,sigma_c) #Minimum symmetrization
  #Omega_hat_US
  sdar_F_sigma3<-matrix_generator3(sigma,sigma_c) #No symmetrization
  TT <- matrixMultiply(t(sdar_F_sigma3),X_plus2)
  #T_hat
  sdar_T<-sdar_F_sigma3+t(sdar_F_sigma3)-matrixMultiply(TT,sdar_F_sigma3) #Desparsified estimator
  
  sdar_F_lowertri1 <- sdar_F_sigma1
  sdar_F_lowertri1[upper.tri(sdar_F_lowertri1)] <- 0
  sdar_F_lowertri1_nodiag <- matrix(0,nrow=p-1,ncol=p)
  sdar_F_lowertri1_nodiag[!upper.tri(sdar_F_lowertri1_nodiag)] <- sdar_F_lowertri1[lower.tri(sdar_F_lowertri1)]
  
  sdar_F_lowertri3 <- sdar_F_sigma3
  sdar_F_lowertri3[upper.tri(sdar_F_lowertri3)] <- 0
  sdar_F_lowertri3_nodiag <- matrix(0,nrow=p-1,ncol=p)
  sdar_F_lowertri3_nodiag[!upper.tri(sdar_F_lowertri3_nodiag)] <- sdar_F_lowertri3[lower.tri(sdar_F_lowertri3)]
  
  #T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S))
  j <- 0
  zvalue_unbias_min <- {}
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
      sigmahat_ij_lowertri <- sqrt((solve(X_plus2[Ajhat_plus,Ajhat_plus])[iAjhat,iAjhat])*(solve(X_plus2[Ajhat_plus,Ajhat_plus])[jAjhat,jAjhat])+(solve(X_plus2[Ajhat_plus,Ajhat_plus])[iAjhat,jAjhat])^2)
      zvalue_ij <- sqrt(n)*theta_hat_ij/sigmahat_ij_lowertri
      zvalue_unbias_min <- c(zvalue_unbias_min,zvalue_ij)
    }
  }
  Ahat_unbias_min <- which(sdar_F_lowertri1_nodiag!=0)
  critical_unbias_min <- Ahat_unbias_min[adaptZ.func(zvalue_unbias_min, alpha)$ac]
  sigma_c_star_unbias_min <- sdar_F_lowertri1_nodiag
  sigma_c_star_unbias_min[critical_unbias_min] <- 0
  sigma_c_star_unbias_min[upper.tri(sigma_c_star_unbias_min)] <- t(sigma_c_star_unbias_min)[!lower.tri(t(sigma_c_star_unbias_min))]
  #T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S))
  sdar_T_star_unbias_min <- matrix_generator3(sigma,sigma_c_star_unbias_min) 
  
  #T(Omega_hat_S|Z0(T_hat),SL(Omega_hat_S)) and T(T_hat|Z0(T_hat),SL(Omega_hat_S))
  zvalue_desparsified_min <- computeZValueDesparsifiedGaussian(sdar_F_lowertri1,sdar_T,sdar_F_sigma3,sigma,n,p)
  Ahat_desparsified_min <- which(sdar_F_lowertri1_nodiag!=0)
  critical_desparsified_min <- Ahat_desparsified_min[adaptZ.func(zvalue_desparsified_min, alpha)$ac]
  sigma_c_star_desparsified_min <- sdar_F_lowertri1_nodiag
  sigma_c_star_desparsified_min[critical_desparsified_min] <- 0
  sigma_c_star_desparsified_min[upper.tri(sigma_c_star_desparsified_min)] <- t(sigma_c_star_desparsified_min)[!lower.tri(t(sigma_c_star_desparsified_min))]
  #T(Omega_hat_S|Z0(T_hat),SL(Omega_hat_S))
  sdar_T_star_desparsified_min_method1 <- matrix_generator3(sigma,sigma_c_star_desparsified_min) 
  Ahat_star_desparsified_min <- which(sdar_T_star_desparsified_min_method1==0)
  #T(T_hat|Z0(T_hat),SL(Omega_hat_S))
  sdar_T_star_desparsified_min_method2 <- sdar_T
  sdar_T_star_desparsified_min_method2[Ahat_star_desparsified_min] <- 0 
  
  #T(T_hat|Z0(T_hat),SL(T_hat))
  zvalue_low <- computeZValueLowGaussian(sdar_T,sdar_F_sigma3,n,p)
  adaptZ <- adaptZ.func(zvalue_low, alpha)
  pvalue_low_adjust <- adaptZ$lfdr
  pmatrix <- matrix(0,nrow=p,ncol=p)
  pmatrix[lower.tri(pmatrix)] <- pvalue_low_adjust
  pmatrix[upper.tri(pmatrix)] <- t(pmatrix)[upper.tri(t(pmatrix))]
  critical_method3 <- which(pmatrix>adaptZ$threshold)
  #T(T_hat|Z0(T_hat),SL(T_hat))
  sdar_T_star_method3 <- sdar_T
  sdar_T_star_method3[critical_method3] <- 0
  
  return(list(Omega_hat_S=sdar_F_sigma1,Omega_hat_star1=sdar_T_star_unbias_min,
              Omega_hat_star2=sdar_T_star_desparsified_min_method1,
              Omega_hat_star3=sdar_T_star_desparsified_min_method2,
              Omega_hat_star4=sdar_T_star_method3,
              T_hat=sdar_T))
}







