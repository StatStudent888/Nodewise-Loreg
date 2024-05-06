seed <- scan("seed.txt",what=integer())
#Sample size
n <- 200 #400
#Dimension
p <- 200 #400

l0time_Gaussian <- {}
l0time_subGaussian <- {}

for(i in seed){
  dir <- paste("L0_time_output/L0_band_n",n,"_p",p,"_seed",i,"_time.RData",sep="")
  load(file=dir)
  
  l0time_Gaussian <- rbind(l0time_Gaussian,timek_Gaussian)
  l0time_subGaussian <- rbind(l0time_subGaussian,timek_subGaussian)
  
  rm(data_Gaussian,data_subGaussian,sdar_F_sigma1_Gaussian,sdar_F_sigma1_subGaussian,timek_Gaussian,timek_subGaussian)
}

#Compute mean of time under 100 replications
l0time_Gaussian_mean <- apply(l0time_Gaussian,2,mean)
l0time_subGaussian_mean <- apply(l0time_subGaussian,2,mean)

#Compute sd of time under 100 replications
l0time_Gaussian_sd <- apply(l0time_Gaussian,2,sd)
l0time_subGaussian_sd <- apply(l0time_subGaussian,2,sd)

l0time_Gaussian_mean
l0time_Gaussian_sd

l0time_subGaussian_mean
l0time_subGaussian_sd

save.image(paste("l0_Gaussian_subGaussian_band_n",n,"_p",p,"_time.RData",sep=""))
