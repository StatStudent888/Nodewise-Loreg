seed <- scan("seed.txt",what=integer())
n <- 200
p <- 200

precision_subGaussian_100 <- {}
sensitivity_subGaussian_100 <- {}
specificity_subGaussian_100 <- {}
MCC_subGaussian_100 <- {}

norm1_subGaussian_100 <- {}
norm2_subGaussian_100 <- {}
normF_subGaussian_100 <- {}
normmax_subGaussian_100 <- {}

MIOtime_subGaussian <- {}

for(i in seed){
  dir <- paste("MIO_subGaussian/MIO_subGaussian_band_n",n,"_p",p,"_seed",i,".RData",sep="")
  load(file=dir)
  
  precision_subGaussian_100 <- c(precision_subGaussian_100,precision)
  sensitivity_subGaussian_100 <- c(sensitivity_subGaussian_100,sensitivity)
  specificity_subGaussian_100 <- c(specificity_subGaussian_100,specificity)
  MCC_subGaussian_100 <- c(MCC_subGaussian_100,MCC)
  
  norm1_subGaussian_100 <- c(norm1_subGaussian_100,MIO_norm1)
  norm2_subGaussian_100 <- c(norm2_subGaussian_100,MIO_norm2)
  normF_subGaussian_100 <- c(normF_subGaussian_100,MIO_normF)
  normmax_subGaussian_100 <- c(normmax_subGaussian_100,MIO_normmax)
  
  MIOtime_subGaussian <- rbind(MIOtime_subGaussian,MIO_time)
  
  rm(precision,sensitivity,specificity,MCC,MIO_norm1,MIO_norm2,MIO_normF,MIO_normmax,MIO_time)
}

precision_subGaussian_mean <- mean(precision_subGaussian_100)
sensitivity_subGaussian_mean <- mean(sensitivity_subGaussian_100)
specificity_subGaussian_mean <- mean(specificity_subGaussian_100)
MCC_subGaussian_mean <- mean(MCC_subGaussian_100)

norm1_subGaussian_mean <- mean(norm1_subGaussian_100)
norm2_subGaussian_mean <- mean(norm2_subGaussian_100)
normF_subGaussian_mean <- mean(normF_subGaussian_100)
normmax_subGaussian_mean <- mean(normmax_subGaussian_100)

MIOtime_subGaussian_mean <- apply(MIOtime_subGaussian,2,mean)

precision_subGaussian_sd <- sd(precision_subGaussian_100)
sensitivity_subGaussian_sd <- sd(sensitivity_subGaussian_100)
specificity_subGaussian_sd <- sd(specificity_subGaussian_100)
MCC_subGaussian_sd <- sd(MCC_subGaussian_100)

norm1_subGaussian_sd <- sd(norm1_subGaussian_100)
norm2_subGaussian_sd <- sd(norm2_subGaussian_100)
normF_subGaussian_sd <- sd(normF_subGaussian_100)
normmax_subGaussian_sd <- sd(normmax_subGaussian_100)

MIOtime_subGaussian_sd <- apply(MIOtime_subGaussian,2,sd)

x1 <- paste(sprintf("%0.3f",precision_subGaussian_mean),"(",sprintf("%0.3f",precision_subGaussian_sd),")"," & ",
            sprintf("%0.3f",sensitivity_subGaussian_mean),"(",sprintf("%0.3f",sensitivity_subGaussian_sd),")"," & ",
            sprintf("%0.3f",specificity_subGaussian_mean),"(",sprintf("%0.3f",specificity_subGaussian_sd),")"," & ",
            sprintf("%0.3f",MCC_subGaussian_mean),"(",sprintf("%0.3f",MCC_subGaussian_sd),")",sep="")

x2 <- paste(sprintf("%0.3f",norm1_subGaussian_mean),"(",sprintf("%0.3f",norm1_subGaussian_sd),")"," & ",
            sprintf("%0.3f",norm2_subGaussian_mean),"(",sprintf("%0.3f",norm2_subGaussian_sd),")"," & ",
            sprintf("%0.3f",normF_subGaussian_mean),"(",sprintf("%0.3f",normF_subGaussian_sd),")"," & ",
            sprintf("%0.3f",normmax_subGaussian_mean),"(",sprintf("%0.3f",normmax_subGaussian_sd),")",sep="")

x2

x1

MIOtime_subGaussian_mean
MIOtime_subGaussian_sd


