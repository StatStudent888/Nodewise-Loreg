seed <- scan("realdata_seed.txt")

p <- 1400

precision_100 <- {}
sensitivity_100 <- {}
specificity_100 <- {}
MCC_100 <- {}
sparsity_100 <- {}
sdar_time_100 <- {}

for(i in seed){
  dir <- paste("realdata_l0/L0_real_data_p",p,"_seed",i,".RData",sep="")
  load(file=dir)
  
  precision_100 <- c(precision_100,precision)
  sensitivity_100 <- c(sensitivity_100,sensitivity)
  specificity_100 <- c(specificity_100,specificity)
  MCC_100 <- c(MCC_100,MCC)
  sparsity_100 <- c(sparsity_100,sparsity)
  sdar_time_100 <- rbind(sdar_time_100,sdar_time)
  
  rm(sdar_time,sparsity,precision,specificity,sensitivity,MCC)
}



precision_mean <- mean(precision_100)
sensitivity_mean <- mean(sensitivity_100)
specificity_mean <- mean(specificity_100)
MCC_mean <- mean(MCC_100)
sparsity_mean <- mean(sparsity_100)
sdar_time_mean <- apply(sdar_time_100,2,mean)

precision_sd <- sd(precision_100)
sensitivity_sd <- sd(sensitivity_100)
specificity_sd <- sd(specificity_100)
MCC_sd <- sd(MCC_100)
sparsity_sd <- sd(sparsity_100)
sdar_time_sd <- apply(sdar_time_100,2,sd)

x1 <- paste(sprintf("%0.3f",precision_mean),"(",sprintf("%0.3f",precision_sd),")"," & ",
            sprintf("%0.3f",sensitivity_mean),"(",sprintf("%0.3f",sensitivity_sd),")"," & ",
            sprintf("%0.3f",specificity_mean),"(",sprintf("%0.3f",specificity_sd),")"," & ",
            sprintf("%0.3f",MCC_mean),"(",sprintf("%0.3f",MCC_sd),")"," & ",
            sprintf("%0.3f",sparsity_mean),"(",sprintf("%0.3f",sparsity_sd),")",sep="")

x1

sdar_time_mean

sdar_time_sd




