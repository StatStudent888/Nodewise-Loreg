seed <- scan("seed.txt",what=integer())
#Sample size
n <- 200 #400
#Dimension
p <- 200 #400

l1time_Gaussian <- {}
l1time_subGaussian <- {}

for(i in seed){
  dir <- paste("L1_time_output/L1_band_n",n,"_p",p,"_seed",i,"_time.RData",sep="")
  load(file=dir)
  
  l1time_Gaussian <- rbind(l1time_Gaussian,timek_Gaussian)
  l1time_subGaussian <- rbind(l1time_subGaussian,timek_subGaussian)
  
  rm(data_Gaussian,data_subGaussian,NS_F_sigma1_Gaussian,NS_F_sigma1_subGaussian,timek_Gaussian,timek_subGaussian)
}

#Compute mean of time under 100 replications
l1time_Gaussian_mean <- apply(l1time_Gaussian,2,mean)
l1time_subGaussian_mean <- apply(l1time_subGaussian,2,mean)

#Compute sd of time under 100 replications
l1time_Gaussian_sd <- apply(l1time_Gaussian,2,sd)
l1time_subGaussian_sd <- apply(l1time_subGaussian,2,sd)

l1time_Gaussian_mean
l1time_Gaussian_sd

l1time_subGaussian_mean
l1time_subGaussian_sd

save.image(paste("l1_Gaussian_subGaussian_band_n",n,"_p",p,"_time.RData",sep=""))