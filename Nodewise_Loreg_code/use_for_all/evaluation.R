#Support recovery performance in simulation
####################################
#theta-------Ture precision matrix
#theta_hat---Estimator of precision matrix
##################################
#A function returns the support recovery performance of estimators on off-diagonal entries
evaluation <- function(theta,theta_hat){
  diag(theta) <- NA
  diag(theta_hat) <- NA
  theta <- matrix(theta[which(!is.na(theta))],nrow=nrow(theta)-1)
  theta_hat <- matrix(theta_hat[which(!is.na(theta_hat))],nrow=nrow(theta_hat)-1)
  TP <- sum((abs(theta)>0.00001)*(theta_hat!=0))
  FP <- sum((abs(theta)<=0.00001)*(theta_hat!=0))
  FN <- sum((abs(theta)>0.00001)*(theta_hat==0))
  TN <- sum((abs(theta)<=0.00001)*(theta_hat==0))
  sensitivity <- TP/(TP+FN)
  specificity <- TN/(FP+TN)
  precision <- TP/(TP+FP)
  MCC <- (TP*TN-FP*FN)/(sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN))
  FNR <- FN/(FN+TN)
  F1 <- 2*precision*sensitivity/(precision+sensitivity)
  return(list(sensitivity=sensitivity,specificity=specificity,MCC=MCC,precision=precision,FNR=FNR,F1=F1))
}
