#Generate symmetric matrix
####################################
#sigma-------A vector contains diagonal elements of a matrix
#sigma_c-----A (p-1)*p matrix that does not contain diagonal elements
##################################
#Minimum symmetrization
matrix_generator1 <- function(sigma,sigma_c)
{
  F_sigma<-diag(sigma) 
  F_sigma[upper.tri(F_sigma,diag=F)]<-sigma_c[upper.tri(sigma_c,diag=F)]
  F_sigma[lower.tri(F_sigma,diag=F)]<-sigma_c[lower.tri(sigma_c,diag=T)]
  F_sigma<-F_sigma 
  F_sigma1<- F_sigma
  F_sigma2<- F_sigma   
  F_sigma1[lower.tri(F_sigma1)] <- t(F_sigma1)[lower.tri(F_sigma1)]
  aF_sigma1<-abs(F_sigma1)
  F_sigma2[upper.tri(F_sigma2)] <- t(F_sigma2)[upper.tri(F_sigma2)]
  aF_sigma2<-abs(F_sigma2)
  aF_sigma<-aF_sigma1-aF_sigma2
  F_sigma[which(aF_sigma<0)]<-F_sigma1[which(aF_sigma<0)]
  F_sigma[which(aF_sigma>0)]<-F_sigma2[which(aF_sigma>0)]
  F_sigma<-as.matrix(F_sigma)
  return(F_sigma)
}

#No symmetrization
matrix_generator3 <- function(sigma,sigma_c){
  F_sigma<-diag(sigma) 
  F_sigma[upper.tri(F_sigma,diag=F)]<-sigma_c[upper.tri(sigma_c,diag=F)]
  F_sigma[lower.tri(F_sigma,diag=F)]<-sigma_c[lower.tri(sigma_c,diag=T)]
  F_sigma<-as.matrix(F_sigma)
  return(F_sigma)
}

