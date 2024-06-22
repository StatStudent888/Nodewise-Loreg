library(huge)
source("use_for_all/matrix_generator.R")
source("use_for_all/omega_generator.R")
source("use_for_all/evaluation.R")
seed <- scan("seed.txt",what=integer())
#Sample size
n=200 #400
#Dimension
p=200 #400

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

precision_Gaussian_100 <- {}
sensitivity_Gaussian_100 <- {}
specificity_Gaussian_100 <- {}
MCC_Gaussian_100 <- {}
precision_subGaussian_100 <- {}
sensitivity_subGaussian_100 <- {}
specificity_subGaussian_100 <- {}
MCC_subGaussian_100 <- {}

norm1_Gaussian_100 <- {}
norm1_subGaussian_100 <- {}
norm2_Gaussian_100 <- {}
norm2_subGaussian_100 <- {}
normF_Gaussian_100 <- {}
normF_subGaussian_100 <- {}
normmax_Gaussian_100 <- {}
normmax_subGaussian_100 <- {}

Glassotime_Gaussian <- {}
Glassotime_subGaussian <- {}


for(i in seed){
  dir <- paste("Glasso_time_output/Glasso_band_n",n,"_p",p,"_seed",i,".RData",sep="")
  load(file=dir)
  
  eva_Gaussian <- evaluation(omega4,Glasso_F_sigma_Gaussian)
  eva_subGaussian <- evaluation(omega4,Glasso_F_sigma_subGaussian)
  precision_Gaussian_100 <- c(precision_Gaussian_100,eva_Gaussian$precision)
  sensitivity_Gaussian_100 <- c(sensitivity_Gaussian_100,eva_Gaussian$sensitivity)
  specificity_Gaussian_100 <- c(specificity_Gaussian_100,eva_Gaussian$specificity)
  MCC_Gaussian_100 <- c(MCC_Gaussian_100,eva_Gaussian$MCC)
  precision_subGaussian_100 <- c(precision_subGaussian_100,eva_subGaussian$precision)
  sensitivity_subGaussian_100 <- c(sensitivity_subGaussian_100,eva_subGaussian$sensitivity)
  specificity_subGaussian_100 <- c(specificity_subGaussian_100,eva_subGaussian$specificity)
  MCC_subGaussian_100 <- c(MCC_subGaussian_100,eva_subGaussian$MCC)
  
  norm1_Gaussian <- norm(Glasso_F_sigma_Gaussian-omega4,"1")
  norm1_subGaussian <- norm(Glasso_F_sigma_subGaussian-omega4,"1")
  norm1_Gaussian_100 <- c(norm1_Gaussian_100,norm1_Gaussian)
  norm1_subGaussian_100 <- c(norm1_subGaussian_100,norm1_subGaussian)
  norm2_Gaussian <- norm(Glasso_F_sigma_Gaussian-omega4,"2")
  norm2_subGaussian <- norm(Glasso_F_sigma_subGaussian-omega4,"2")
  norm2_Gaussian_100 <- c(norm2_Gaussian_100,norm2_Gaussian)
  norm2_subGaussian_100 <- c(norm2_subGaussian_100,norm2_subGaussian)
  normF_Gaussian <- norm(Glasso_F_sigma_Gaussian-omega4,"F")
  normF_subGaussian <- norm(Glasso_F_sigma_subGaussian-omega4,"F")
  normF_Gaussian_100 <- c(normF_Gaussian_100,normF_Gaussian)
  normF_subGaussian_100 <- c(normF_subGaussian_100,normF_subGaussian)
  normmax_Gaussian <- norm(Glasso_F_sigma_Gaussian-omega4,"M")
  normmax_subGaussian <- norm(Glasso_F_sigma_subGaussian-omega4,"M")
  normmax_Gaussian_100 <- c(normmax_Gaussian_100,normmax_Gaussian)
  normmax_subGaussian_100 <- c(normmax_subGaussian_100,normmax_subGaussian)
  
  Glassotime_Gaussian <- rbind(Glassotime_Gaussian,timek_Gaussian)
  Glassotime_subGaussian <- rbind(Glassotime_subGaussian,timek_subGaussian)
  
  rm(data_Gaussian,data_subGaussian,Glasso_F_sigma_Gaussian,Glasso_F_sigma_subGaussian,timek_Gaussian,timek_subGaussian)
}


precision_Gaussian_mean <- mean(precision_Gaussian_100)
sensitivity_Gaussian_mean <- mean(sensitivity_Gaussian_100)
specificity_Gaussian_mean <- mean(specificity_Gaussian_100)
MCC_Gaussian_mean <- mean(MCC_Gaussian_100)
precision_subGaussian_mean <- mean(precision_subGaussian_100)
sensitivity_subGaussian_mean <- mean(sensitivity_subGaussian_100)
specificity_subGaussian_mean <- mean(specificity_subGaussian_100)
MCC_subGaussian_mean <- mean(MCC_subGaussian_100)

norm1_Gaussian_mean <- mean(norm1_Gaussian_100)
norm1_subGaussian_mean <- mean(norm1_subGaussian_100)
norm2_Gaussian_mean <- mean(norm2_Gaussian_100)
norm2_subGaussian_mean <- mean(norm2_subGaussian_100)
normF_Gaussian_mean <- mean(normF_Gaussian_100)
normF_subGaussian_mean <- mean(normF_subGaussian_100)
normmax_Gaussian_mean <- mean(normmax_Gaussian_100)
normmax_subGaussian_mean <- mean(normmax_subGaussian_100)

Glassotime_Gaussian_mean <- apply(Glassotime_Gaussian,2,mean)
Glassotime_subGaussian_mean <- apply(Glassotime_subGaussian,2,mean)

precision_Gaussian_sd <- sd(precision_Gaussian_100)
sensitivity_Gaussian_sd <- sd(sensitivity_Gaussian_100)
specificity_Gaussian_sd <- sd(specificity_Gaussian_100)
MCC_Gaussian_sd <- sd(MCC_Gaussian_100)
precision_subGaussian_sd <- sd(precision_subGaussian_100)
sensitivity_subGaussian_sd <- sd(sensitivity_subGaussian_100)
specificity_subGaussian_sd <- sd(specificity_subGaussian_100)
MCC_subGaussian_sd <- sd(MCC_subGaussian_100)

norm1_Gaussian_sd <- sd(norm1_Gaussian_100)
norm1_subGaussian_sd <- sd(norm1_subGaussian_100)
norm2_Gaussian_sd <- sd(norm2_Gaussian_100)
norm2_subGaussian_sd <- sd(norm2_subGaussian_100)
normF_Gaussian_sd <- sd(normF_Gaussian_100)
normF_subGaussian_sd <- sd(normF_subGaussian_100)
normmax_Gaussian_sd <- sd(normmax_Gaussian_100)
normmax_subGaussian_sd <- sd(normmax_subGaussian_100)

Glassotime_Gaussian_sd <- apply(Glassotime_Gaussian,2,sd)
Glassotime_subGaussian_sd <- apply(Glassotime_subGaussian,2,sd)

x1 <- paste(sprintf("%0.3f",precision_Gaussian_mean),"(",sprintf("%0.3f",precision_Gaussian_sd),")"," & ",
            sprintf("%0.3f",sensitivity_Gaussian_mean),"(",sprintf("%0.3f",sensitivity_Gaussian_sd),")"," & ",
            sprintf("%0.3f",specificity_Gaussian_mean),"(",sprintf("%0.3f",specificity_Gaussian_sd),")"," & ",
            sprintf("%0.3f",MCC_Gaussian_mean),"(",sprintf("%0.3f",MCC_Gaussian_sd),")",sep="")

x2 <- paste(sprintf("%0.3f",norm1_Gaussian_mean),"(",sprintf("%0.3f",norm1_Gaussian_sd),")"," & ",
            sprintf("%0.3f",norm2_Gaussian_mean),"(",sprintf("%0.3f",norm2_Gaussian_sd),")"," & ",
            sprintf("%0.3f",normF_Gaussian_mean),"(",sprintf("%0.3f",normF_Gaussian_sd),")"," & ",
            sprintf("%0.3f",normmax_Gaussian_mean),"(",sprintf("%0.3f",normmax_Gaussian_sd),")",sep="")

x3 <- paste(sprintf("%0.3f",precision_subGaussian_mean),"(",sprintf("%0.3f",precision_subGaussian_sd),")"," & ",
            sprintf("%0.3f",sensitivity_subGaussian_mean),"(",sprintf("%0.3f",sensitivity_subGaussian_sd),")"," & ",
            sprintf("%0.3f",specificity_subGaussian_mean),"(",sprintf("%0.3f",specificity_subGaussian_sd),")"," & ",
            sprintf("%0.3f",MCC_subGaussian_mean),"(",sprintf("%0.3f",MCC_subGaussian_sd),")",sep="")

x4 <- paste(sprintf("%0.3f",norm1_subGaussian_mean),"(",sprintf("%0.3f",norm1_subGaussian_sd),")"," & ",
            sprintf("%0.3f",norm2_subGaussian_mean),"(",sprintf("%0.3f",norm2_subGaussian_sd),")"," & ",
            sprintf("%0.3f",normF_subGaussian_mean),"(",sprintf("%0.3f",normF_subGaussian_sd),")"," & ",
            sprintf("%0.3f",normmax_subGaussian_mean),"(",sprintf("%0.3f",normmax_subGaussian_sd),")",sep="")

x1

x2

Glassotime_Gaussian_mean
Glassotime_Gaussian_sd

x3

x4

Glassotime_subGaussian_mean
Glassotime_subGaussian_sd

save.image(paste("Glasso_Gaussian_subGaussian_band_n",n,"_p",p,"_all.RData",sep=""))