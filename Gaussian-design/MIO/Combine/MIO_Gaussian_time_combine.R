seed <- scan("seed.txt",what=integer())
n <- 200
p <- 200

precision_Gaussian_100 <- {}
sensitivity_Gaussian_100 <- {}
specificity_Gaussian_100 <- {}
MCC_Gaussian_100 <- {}

norm1_Gaussian_100 <- {}
norm2_Gaussian_100 <- {}
normF_Gaussian_100 <- {}
normmax_Gaussian_100 <- {}

MIOtime_Gaussian <- {}

for(i in seed){
  dir <- paste("MIO_Gaussian/MIO_Gaussian_band_n",n,"_p",p,"_seed",i,".RData",sep="")
  load(file=dir)
  
  precision_Gaussian_100 <- c(precision_Gaussian_100,precision)
  sensitivity_Gaussian_100 <- c(sensitivity_Gaussian_100,sensitivity)
  specificity_Gaussian_100 <- c(specificity_Gaussian_100,specificity)
  MCC_Gaussian_100 <- c(MCC_Gaussian_100,MCC)
  
  norm1_Gaussian_100 <- c(norm1_Gaussian_100,MIO_norm1)
  norm2_Gaussian_100 <- c(norm2_Gaussian_100,MIO_norm2)
  normF_Gaussian_100 <- c(normF_Gaussian_100,MIO_normF)
  normmax_Gaussian_100 <- c(normmax_Gaussian_100,MIO_normmax)
  
  MIOtime_Gaussian <- rbind(MIOtime_Gaussian,MIO_time)
  
  rm(precision,sensitivity,specificity,MCC,MIO_norm1,MIO_norm2,MIO_normF,MIO_normmax,MIO_time)
}


precision_Gaussian_mean <- mean(precision_Gaussian_100)
sensitivity_Gaussian_mean <- mean(sensitivity_Gaussian_100)
specificity_Gaussian_mean <- mean(specificity_Gaussian_100)
MCC_Gaussian_mean <- mean(MCC_Gaussian_100)

norm1_Gaussian_mean <- mean(norm1_Gaussian_100)
norm2_Gaussian_mean <- mean(norm2_Gaussian_100)
normF_Gaussian_mean <- mean(normF_Gaussian_100)
normmax_Gaussian_mean <- mean(normmax_Gaussian_100)

MIOtime_Gaussian_mean <- apply(MIOtime_Gaussian,2,mean)

precision_Gaussian_sd <- sd(precision_Gaussian_100)
sensitivity_Gaussian_sd <- sd(sensitivity_Gaussian_100)
specificity_Gaussian_sd <- sd(specificity_Gaussian_100)
MCC_Gaussian_sd <- sd(MCC_Gaussian_100)

norm1_Gaussian_sd <- sd(norm1_Gaussian_100)
norm2_Gaussian_sd <- sd(norm2_Gaussian_100)
normF_Gaussian_sd <- sd(normF_Gaussian_100)
normmax_Gaussian_sd <- sd(normmax_Gaussian_100)

MIOtime_Gaussian_sd <- apply(MIOtime_Gaussian,2,sd)

x1 <- paste(sprintf("%0.3f",precision_Gaussian_mean),"(",sprintf("%0.3f",precision_Gaussian_sd),")"," & ",
            sprintf("%0.3f",sensitivity_Gaussian_mean),"(",sprintf("%0.3f",sensitivity_Gaussian_sd),")"," & ",
            sprintf("%0.3f",specificity_Gaussian_mean),"(",sprintf("%0.3f",specificity_Gaussian_sd),")"," & ",
            sprintf("%0.3f",MCC_Gaussian_mean),"(",sprintf("%0.3f",MCC_Gaussian_sd),")",sep="")

x2 <- paste(sprintf("%0.3f",norm1_Gaussian_mean),"(",sprintf("%0.3f",norm1_Gaussian_sd),")"," & ",
            sprintf("%0.3f",norm2_Gaussian_mean),"(",sprintf("%0.3f",norm2_Gaussian_sd),")"," & ",
            sprintf("%0.3f",normF_Gaussian_mean),"(",sprintf("%0.3f",normF_Gaussian_sd),")"," & ",
            sprintf("%0.3f",normmax_Gaussian_mean),"(",sprintf("%0.3f",normmax_Gaussian_sd),")",sep="")

x2

x1

MIOtime_Gaussian_mean
MIOtime_Gaussian_sd

