####################################
#x-----------Data of a given testing subject
#sigma-------Estimator of precision matrix
#miu1--------Sample mean of the first class
#miu2--------Sample mean of the second class
#pi1---------Proportion of the first class
#pi2---------Proportion of the second class
##################################
#A function returns LDA score
score <- function(x,sigma,miu1,miu2,pi1,pi2){
  derta1 <- t(x)%*%sigma%*%miu1-0.5*t(miu1)%*%sigma%*%miu1+log(pi1)
  derta2 <- t(x)%*%sigma%*%miu2-0.5*t(miu2)%*%sigma%*%miu2+log(pi2)
  return(list(derta1=derta1,derta2=derta2))
}
