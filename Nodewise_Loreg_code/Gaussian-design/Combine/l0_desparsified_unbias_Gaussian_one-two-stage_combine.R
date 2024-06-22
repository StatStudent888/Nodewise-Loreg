library(huge)
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
seed <- scan("seed.txt",what=integer())
#Sample size
n <- 200 #400
#Dimension
p <- 200 #400

#Generate Band graph
set.seed(-1)
omega2 <- omega_generator1(p,1,0.5,0.3)
Sigma2 <- solve(omega2)

#Generate Random graph
#graph= "random"
#set.seed(-1)
#generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, prob=4/p)
#omega2 <- generator$omega
#Sigma2 <- generator$sigma

#Generate Hub graph
#graph= "hub"
#set.seed(-1)
#generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
#omega2 <- generator$omega
#Sigma2 <- generator$sigma

#Generate Cluster graph
#graph= "cluster"
#set.seed(-1)
#generator<-huge.generator(n, d=p, graph, verbose=FALSE, v=1, u=0, g=p/10)
#omega2 <- generator$omega
#Sigma2 <- generator$sigma

precision_min_100 <- {}
sensitivity_min_100 <- {}
specificity_min_100 <- {}
MCC_min_100 <- {}
precision_star_desparsified_min_100 <- {}
sensitivity_star_desparsified_min_100 <- {}
specificity_star_desparsified_min_100 <- {}
MCC_star_desparsified_min_100 <- {}
precision_star_unbias_min_100 <- {}
sensitivity_star_unbias_min_100 <- {}
specificity_star_unbias_min_100 <- {}
MCC_star_unbias_min_100 <- {}

precision_star_method3_100 <- {}
sensitivity_star_method3_100 <- {}
specificity_star_method3_100 <- {}
MCC_star_method3_100 <- {}

norm1_min_100 <- {}
norm2_min_100 <- {}
normF_min_100 <- {}
normmax_min_100 <- {}
norm1_star_desparsified_min_method1_100 <- {}
norm2_star_desparsified_min_method1_100 <- {}
normF_star_desparsified_min_method1_100 <- {}
normmax_star_desparsified_min_method1_100 <- {}
norm1_star_unbias_min_100 <- {}
norm2_star_unbias_min_100 <- {}
normF_star_unbias_min_100 <- {}
normmax_star_unbias_min_100 <- {}

norm1_star_desparsified_min_method2_100 <- {}
norm2_star_desparsified_min_method2_100 <- {}
normF_star_desparsified_min_method2_100 <- {}
normmax_star_desparsified_min_method2_100 <- {}

norm1_star_method3_100 <- {}
norm2_star_method3_100 <- {}
normF_star_method3_100 <- {}
normmax_star_method3_100 <- {}

norm1_T_100 <- {}
norm2_T_100 <- {}
normF_T_100 <- {}
normmax_T_100 <- {}

for(i in seed){
  dir <- paste("l0_desparsified_unbias_Gaussian_one-two-stage_lfdr_output/L0_desparsified_unbias_Gaussian_band_n",n,"_p",p,"_seed",i,"_one-two-stage_lfdr.RData",sep="")
  load(file=dir)
  
  precision_min_100 <- c(precision_min_100,precision_min)
  sensitivity_min_100 <- c(sensitivity_min_100,sensitivity_min)
  specificity_min_100 <- c(specificity_min_100,specificity_min)
  MCC_min_100 <- c(MCC_min_100,MCC_min)
  precision_star_desparsified_min_100 <- c(precision_star_desparsified_min_100,precision_star_desparsified_min)
  sensitivity_star_desparsified_min_100 <- c(sensitivity_star_desparsified_min_100,sensitivity_star_desparsified_min)
  specificity_star_desparsified_min_100 <- c(specificity_star_desparsified_min_100,specificity_star_desparsified_min)
  MCC_star_desparsified_min_100 <- c(MCC_star_desparsified_min_100,MCC_star_desparsified_min)
  precision_star_unbias_min_100 <- c(precision_star_unbias_min_100,precision_star_unbias_min)
  sensitivity_star_unbias_min_100 <- c(sensitivity_star_unbias_min_100,sensitivity_star_unbias_min)
  specificity_star_unbias_min_100 <- c(specificity_star_unbias_min_100,specificity_star_unbias_min)
  MCC_star_unbias_min_100 <- c(MCC_star_unbias_min_100,MCC_star_unbias_min)
  
  precision_star_method3_100 <- c(precision_star_method3_100,precision_star_method3)
  sensitivity_star_method3_100 <- c(sensitivity_star_method3_100,sensitivity_star_method3)
  specificity_star_method3_100 <- c(specificity_star_method3_100,specificity_star_method3)
  MCC_star_method3_100 <- c(MCC_star_method3_100,MCC_star_method3)
  
  #T_hat
  norm1_T <- norm(sdar_T-omega2,"1")
  norm2_T <- norm(sdar_T-omega2,"2")
  normF_T <- norm(sdar_T-omega2,"F")
  normmax_T <- norm(sdar_T-omega2,"M")
  norm1_T_100 <- c(norm1_T_100,norm1_T)
  norm2_T_100 <- c(norm2_T_100,norm2_T)
  normF_T_100 <- c(normF_T_100,normF_T)
  normmax_T_100 <- c(normmax_T_100,normmax_T)
  
  norm1_min_100 <- c(norm1_min_100,sdar_min_norm1)
  norm2_min_100 <- c(norm2_min_100,sdar_min_norm2)
  normF_min_100 <- c(normF_min_100,sdar_min_normF)
  normmax_min_100 <- c(normmax_min_100,sdar_min_normmax)
  norm1_star_desparsified_min_method1_100 <- c(norm1_star_desparsified_min_method1_100,sdar_norm1_T_star_desparsified_min_method1)
  norm2_star_desparsified_min_method1_100 <- c(norm2_star_desparsified_min_method1_100,sdar_norm2_T_star_desparsified_min_method1)
  normF_star_desparsified_min_method1_100 <- c(normF_star_desparsified_min_method1_100,sdar_normF_T_star_desparsified_min_method1)
  normmax_star_desparsified_min_method1_100 <- c(normmax_star_desparsified_min_method1_100,sdar_normmax_T_star_desparsified_min_method1)
  norm1_star_unbias_min_100 <- c(norm1_star_unbias_min_100,sdar_norm1_T_star_unbias_min)
  norm2_star_unbias_min_100 <- c(norm2_star_unbias_min_100,sdar_norm2_T_star_unbias_min)
  normF_star_unbias_min_100 <- c(normF_star_unbias_min_100,sdar_normF_T_star_unbias_min)
  normmax_star_unbias_min_100 <- c(normmax_star_unbias_min_100,sdar_normmax_T_star_unbias_min)
  
  norm1_star_desparsified_min_method2_100 <- c(norm1_star_desparsified_min_method2_100,sdar_norm1_T_star_desparsified_min_method2)
  norm2_star_desparsified_min_method2_100 <- c(norm2_star_desparsified_min_method2_100,sdar_norm2_T_star_desparsified_min_method2)
  normF_star_desparsified_min_method2_100 <- c(normF_star_desparsified_min_method2_100,sdar_normF_T_star_desparsified_min_method2)
  normmax_star_desparsified_min_method2_100 <- c(normmax_star_desparsified_min_method2_100,sdar_normmax_T_star_desparsified_min_method2)
  
  norm1_star_method3_100 <- c(norm1_star_method3_100,sdar_norm1_T_star_method3)
  norm2_star_method3_100 <- c(norm2_star_method3_100,sdar_norm2_T_star_method3)
  normF_star_method3_100 <- c(normF_star_method3_100,sdar_normF_T_star_method3)
  normmax_star_method3_100 <- c(normmax_star_method3_100,sdar_normmax_T_star_method3)
  
  rm(data,sdar_F_sigma1,sdar_F_sigma3,sdar_T,
     sdar_T_star_unbias_min,
     sdar_T_star_desparsified_min_method1,sdar_T_star_desparsified_min_method2,
     sdar_T_star_method3,precision_min,precision_star_unbias_min,
     sensitivity_min,sensitivity_star_unbias_min,specificity_min,
     specificity_star_unbias_min,MCC_min,MCC_star_unbias_min,F1_min,F1_star_unbias_min,
     precision_star_desparsified_min,
     sensitivity_star_desparsified_min,specificity_star_desparsified_min,
     MCC_star_desparsified_min,F1_star_desparsified_min,precision_star_method3,
     sensitivity_star_method3,specificity_star_method3,MCC_star_method3,F1_star_method3,
     sdar_min_norm1,sdar_min_norm2,sdar_min_normF,sdar_min_normmax,
     sdar_norm1_T_star_unbias_min,sdar_norm2_T_star_unbias_min,
     sdar_normF_T_star_unbias_min,sdar_normmax_T_star_unbias_min,
     sdar_norm1_T_star_desparsified_min_method1,
     sdar_norm2_T_star_desparsified_min_method1,
     sdar_normF_T_star_desparsified_min_method1,
     sdar_normmax_T_star_desparsified_min_method1,
     sdar_norm1_T_star_desparsified_min_method2,
     sdar_norm2_T_star_desparsified_min_method2,
     sdar_normF_T_star_desparsified_min_method2,
     sdar_normmax_T_star_desparsified_min_method2,
     sdar_norm1_T_star_method3,sdar_norm2_T_star_method3,
     sdar_normF_T_star_method3,sdar_normmax_T_star_method3)
}

#Compute mean of support recovery performance under 100 replications
precision_min_mean <- mean(precision_min_100)
sensitivity_min_mean <- mean(sensitivity_min_100)
specificity_min_mean <- mean(specificity_min_100)
MCC_min_mean <- mean(MCC_min_100)
precision_star_desparsified_min_mean <- mean(precision_star_desparsified_min_100)
sensitivity_star_desparsified_min_mean <- mean(sensitivity_star_desparsified_min_100)
specificity_star_desparsified_min_mean <- mean(specificity_star_desparsified_min_100)
MCC_star_desparsified_min_mean <- mean(MCC_star_desparsified_min_100)
precision_star_unbias_min_mean <- mean(precision_star_unbias_min_100)
sensitivity_star_unbias_min_mean <- mean(sensitivity_star_unbias_min_100)
specificity_star_unbias_min_mean <- mean(specificity_star_unbias_min_100)
MCC_star_unbias_min_mean <- mean(MCC_star_unbias_min_100)

precision_star_method3_mean <- mean(precision_star_method3_100)
sensitivity_star_method3_mean <- mean(sensitivity_star_method3_100)
specificity_star_method3_mean <- mean(specificity_star_method3_100)
MCC_star_method3_mean <- mean(MCC_star_method3_100)

#Compute mean of matrix norm losses under 100 replications
norm1_T_mean <- mean(norm1_T_100)
norm2_T_mean <- mean(norm2_T_100)
normF_T_mean <- mean(normF_T_100)
normmax_T_mean <- mean(normmax_T_100)

norm1_min_mean <- mean(norm1_min_100)
norm2_min_mean <- mean(norm2_min_100)
normF_min_mean <- mean(normF_min_100)
normmax_min_mean <- mean(normmax_min_100)
norm1_star_desparsified_min_method1_mean <- mean(norm1_star_desparsified_min_method1_100)
norm2_star_desparsified_min_method1_mean <- mean(norm2_star_desparsified_min_method1_100)
normF_star_desparsified_min_method1_mean <- mean(normF_star_desparsified_min_method1_100)
normmax_star_desparsified_min_method1_mean <- mean(normmax_star_desparsified_min_method1_100)
norm1_star_unbias_min_mean <- mean(norm1_star_unbias_min_100)
norm2_star_unbias_min_mean <- mean(norm2_star_unbias_min_100)
normF_star_unbias_min_mean <- mean(normF_star_unbias_min_100)
normmax_star_unbias_min_mean <- mean(normmax_star_unbias_min_100)

norm1_star_desparsified_min_method2_mean <- mean(norm1_star_desparsified_min_method2_100)
norm2_star_desparsified_min_method2_mean <- mean(norm2_star_desparsified_min_method2_100)
normF_star_desparsified_min_method2_mean <- mean(normF_star_desparsified_min_method2_100)
normmax_star_desparsified_min_method2_mean <- mean(normmax_star_desparsified_min_method2_100)

norm1_star_method3_mean <- mean(norm1_star_method3_100)
norm2_star_method3_mean <- mean(norm2_star_method3_100)
normF_star_method3_mean <- mean(normF_star_method3_100)
normmax_star_method3_mean <- mean(normmax_star_method3_100)

#Compute sd of support recovery performance under 100 replications
precision_min_sd <- sd(precision_min_100)
sensitivity_min_sd <- sd(sensitivity_min_100)
specificity_min_sd <- sd(specificity_min_100)
MCC_min_sd <- sd(MCC_min_100)
precision_star_desparsified_min_sd <- sd(precision_star_desparsified_min_100)
sensitivity_star_desparsified_min_sd <- sd(sensitivity_star_desparsified_min_100)
specificity_star_desparsified_min_sd <- sd(specificity_star_desparsified_min_100)
MCC_star_desparsified_min_sd <- sd(MCC_star_desparsified_min_100)
precision_star_unbias_min_sd <- sd(precision_star_unbias_min_100)
sensitivity_star_unbias_min_sd <- sd(sensitivity_star_unbias_min_100)
specificity_star_unbias_min_sd <- sd(specificity_star_unbias_min_100)
MCC_star_unbias_min_sd <- sd(MCC_star_unbias_min_100)

precision_star_method3_sd <- sd(precision_star_method3_100)
sensitivity_star_method3_sd <- sd(sensitivity_star_method3_100)
specificity_star_method3_sd <- sd(specificity_star_method3_100)
MCC_star_method3_sd <- sd(MCC_star_method3_100)

#Compute sd of matrix norm losses under 100 replications
norm1_T_sd <- sd(norm1_T_100)
norm2_T_sd <- sd(norm2_T_100)
normF_T_sd <- sd(normF_T_100)
normmax_T_sd <- sd(normmax_T_100)

norm1_min_sd <- sd(norm1_min_100)
norm2_min_sd <- sd(norm2_min_100)
normF_min_sd <- sd(normF_min_100)
normmax_min_sd <- sd(normmax_min_100)
norm1_star_desparsified_min_method1_sd <- sd(norm1_star_desparsified_min_method1_100)
norm2_star_desparsified_min_method1_sd <- sd(norm2_star_desparsified_min_method1_100)
normF_star_desparsified_min_method1_sd <- sd(normF_star_desparsified_min_method1_100)
normmax_star_desparsified_min_method1_sd <- sd(normmax_star_desparsified_min_method1_100)
norm1_star_unbias_min_sd <- sd(norm1_star_unbias_min_100)
norm2_star_unbias_min_sd <- sd(norm2_star_unbias_min_100)
normF_star_unbias_min_sd <- sd(normF_star_unbias_min_100)
normmax_star_unbias_min_sd <- sd(normmax_star_unbias_min_100)

norm1_star_desparsified_min_method2_sd <- sd(norm1_star_desparsified_min_method2_100)
norm2_star_desparsified_min_method2_sd <- sd(norm2_star_desparsified_min_method2_100)
normF_star_desparsified_min_method2_sd <- sd(normF_star_desparsified_min_method2_100)
normmax_star_desparsified_min_method2_sd <- sd(normmax_star_desparsified_min_method2_100)

norm1_star_method3_sd <- sd(norm1_star_method3_100)
norm2_star_method3_sd <- sd(norm2_star_method3_100)
normF_star_method3_sd <- sd(normF_star_method3_100)
normmax_star_method3_sd <- sd(normmax_star_method3_100)


x1 <- paste(sprintf("%0.3f",precision_min_mean),"(",sprintf("%0.3f",precision_min_sd),")"," / ",sprintf("%0.3f",precision_star_unbias_min_mean),"(",sprintf("%0.3f",precision_star_unbias_min_sd),")"," & ",
            sprintf("%0.3f",sensitivity_min_mean),"(",sprintf("%0.3f",sensitivity_min_sd),")"," / ",sprintf("%0.3f",sensitivity_star_unbias_min_mean),"(",sprintf("%0.3f",sensitivity_star_unbias_min_sd),")"," & ",
            sprintf("%0.3f",specificity_min_mean),"(",sprintf("%0.3f",specificity_min_sd),")"," / ",sprintf("%0.3f",specificity_star_unbias_min_mean),"(",sprintf("%0.3f",specificity_star_unbias_min_sd),")"," & ",
            sprintf("%0.3f",MCC_min_mean),"(",sprintf("%0.3f",MCC_min_sd),")"," / ",sprintf("%0.3f",MCC_star_unbias_min_mean),"(",sprintf("%0.3f",MCC_star_unbias_min_sd),")",sep="")
x2 <- paste(sprintf("%0.3f",precision_min_mean),"(",sprintf("%0.3f",precision_min_sd),")"," / ",sprintf("%0.3f",precision_star_desparsified_min_mean),"(",sprintf("%0.3f",precision_star_desparsified_min_sd),")"," & ",
            sprintf("%0.3f",sensitivity_min_mean),"(",sprintf("%0.3f",sensitivity_min_sd),")"," / ",sprintf("%0.3f",sensitivity_star_desparsified_min_mean),"(",sprintf("%0.3f",sensitivity_star_desparsified_min_sd),")"," & ",
            sprintf("%0.3f",specificity_min_mean),"(",sprintf("%0.3f",specificity_min_sd),")"," / ",sprintf("%0.3f",specificity_star_desparsified_min_mean),"(",sprintf("%0.3f",specificity_star_desparsified_min_sd),")"," & ",
            sprintf("%0.3f",MCC_min_mean),"(",sprintf("%0.3f",MCC_min_sd),")"," / ",sprintf("%0.3f",MCC_star_desparsified_min_mean),"(",sprintf("%0.3f",MCC_star_desparsified_min_sd),")",sep="")

x3 <- paste(sprintf("%0.3f",precision_star_method3_mean),"(",sprintf("%0.3f",precision_star_method3_sd),")"," & ",
            sprintf("%0.3f",sensitivity_star_method3_mean),"(",sprintf("%0.3f",sensitivity_star_method3_sd),")"," & ",
            sprintf("%0.3f",specificity_star_method3_mean),"(",sprintf("%0.3f",specificity_star_method3_sd),")"," & ",
            sprintf("%0.3f",MCC_star_method3_mean),"(",sprintf("%0.3f",MCC_star_method3_sd),")",sep="")

x4 <- paste(sprintf("%0.3f",norm1_min_mean),"(",sprintf("%0.3f",norm1_min_sd),")"," / ",sprintf("%0.3f",norm1_star_unbias_min_mean),"(",sprintf("%0.3f",norm1_star_unbias_min_sd),")"," & ",
            sprintf("%0.3f",norm2_min_mean),"(",sprintf("%0.3f",norm2_min_sd),")"," / ",sprintf("%0.3f",norm2_star_unbias_min_mean),"(",sprintf("%0.3f",norm2_star_unbias_min_sd),")"," & ",
            sprintf("%0.3f",normF_min_mean),"(",sprintf("%0.3f",normF_min_sd),")"," / ",sprintf("%0.3f",normF_star_unbias_min_mean),"(",sprintf("%0.3f",normF_star_unbias_min_sd),")"," & ",
            sprintf("%0.3f",normmax_min_mean),"(",sprintf("%0.3f",normmax_min_sd),")"," / ",sprintf("%0.3f",normmax_star_unbias_min_mean),"(",sprintf("%0.3f",normmax_star_unbias_min_sd),")",sep="")
x5 <- paste(sprintf("%0.3f",norm1_min_mean),"(",sprintf("%0.3f",norm1_min_sd),")"," / ",sprintf("%0.3f",norm1_star_desparsified_min_method1_mean),"(",sprintf("%0.3f",norm1_star_desparsified_min_method1_sd),")"," & ",
            sprintf("%0.3f",norm2_min_mean),"(",sprintf("%0.3f",norm2_min_sd),")"," / ",sprintf("%0.3f",norm2_star_desparsified_min_method1_mean),"(",sprintf("%0.3f",norm2_star_desparsified_min_method1_sd),")"," & ",
            sprintf("%0.3f",normF_min_mean),"(",sprintf("%0.3f",normF_min_sd),")"," / ",sprintf("%0.3f",normF_star_desparsified_min_method1_mean),"(",sprintf("%0.3f",normF_star_desparsified_min_method1_sd),")"," & ",
            sprintf("%0.3f",normmax_min_mean),"(",sprintf("%0.3f",normmax_min_sd),")"," / ",sprintf("%0.3f",normmax_star_desparsified_min_method1_mean),"(",sprintf("%0.3f",normmax_star_desparsified_min_method1_sd),")",sep="")

x6 <- paste(sprintf("%0.3f",norm1_star_desparsified_min_method2_mean),"(",sprintf("%0.3f",norm1_star_desparsified_min_method2_sd),")"," & ",
             sprintf("%0.3f",norm2_star_desparsified_min_method2_mean),"(",sprintf("%0.3f",norm2_star_desparsified_min_method2_sd),")"," & ",
             sprintf("%0.3f",normF_star_desparsified_min_method2_mean),"(",sprintf("%0.3f",normF_star_desparsified_min_method2_sd),")"," & ",
             sprintf("%0.3f",normmax_star_desparsified_min_method2_mean),"(",sprintf("%0.3f",normmax_star_desparsified_min_method2_sd),")",sep="")

x7 <- paste(sprintf("%0.3f",norm1_star_method3_mean),"(",sprintf("%0.3f",norm1_star_method3_sd),")"," & ",
             sprintf("%0.3f",norm2_star_method3_mean),"(",sprintf("%0.3f",norm2_star_method3_sd),")"," & ",
             sprintf("%0.3f",normF_star_method3_mean),"(",sprintf("%0.3f",normF_star_method3_sd),")"," & ",
             sprintf("%0.3f",normmax_star_method3_mean),"(",sprintf("%0.3f",normmax_star_method3_sd),")",sep="")

x8 <- paste(sprintf("%0.3f",norm1_T_mean),"(",sprintf("%0.3f",norm1_T_sd),")"," & ",
             sprintf("%0.3f",norm2_T_mean),"(",sprintf("%0.3f",norm2_T_sd),")"," & ",
             sprintf("%0.3f",normF_T_mean),"(",sprintf("%0.3f",normF_T_sd),")"," & ",
             sprintf("%0.3f",normmax_T_mean),"(",sprintf("%0.3f",normmax_T_sd),")",sep="")

x1
x2
x3
x4
x5
x6
x7
x8

save.image(paste("L0_desparsified_unbias_Gaussian_band_n",n,"_p",p,"_one-two-stage_lfdr.RData",sep=""))
