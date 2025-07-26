seed <- scan("seed.txt",what=integer())
n <- 200 #400
p <- 200 #400,1000,4000

l0time_Gaussian <- {}
l0time_subGaussian <- {}

for(i in seed){
  dir <- paste("L0_time_Eigen/L0_band_n",n,"_p",p,"_seed",i,"_time.RData",sep="")
  load(file=dir)
  
  l0time_Gaussian <- rbind(l0time_Gaussian,timek_Gaussian)
  l0time_subGaussian <- rbind(l0time_subGaussian,timek_subGaussian)
  
  rm(timek_Gaussian,timek_subGaussian)
}

l0time_Gaussian_mean <- apply(l0time_Gaussian,2,mean)
l0time_subGaussian_mean <- apply(l0time_subGaussian,2,mean)

l0time_Gaussian_sd <- apply(l0time_Gaussian,2,sd)
l0time_subGaussian_sd <- apply(l0time_subGaussian,2,sd)

l0time_Gaussian_mean
l0time_Gaussian_sd

l0time_subGaussian_mean
l0time_subGaussian_sd


