#Classification performance in real data
####################################
#derta_pre---Estimation of class
#dera-ture---Ture class
##################################
#A function returns classification performance 
real_evaluation <- function(derta_pre,derta_true){
  TP <- as.numeric(sum((derta_pre==1)*(derta_true==1)))
  FP <- as.numeric(sum((derta_pre==1)*(derta_true==0)))
  FN <- as.numeric(sum((derta_pre==0)*(derta_true==1)))
  TN <- as.numeric(sum((derta_pre==0)*(derta_true==0)))
  sensitivity <- TP/(TP+FN)
  specificity <- TN/(FP+TN)
  precision <- TP/(TP+FP)
  MCC <- (TP*TN-FP*FN)/(sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN))
  return(list(sensitivity=sensitivity,specificity=specificity,MCC=MCC,precision=precision))
}


