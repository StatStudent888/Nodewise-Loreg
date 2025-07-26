seed <- scan("seed.txt",what=integer())
n <- 200 #400
p <- 200 #400,1000,4000

l1time_Gaussian <- {}
l1time_subGaussian <- {}

for(i in seed){
  dir <- paste("L1_time/L1_band_n",n,"_p",p,"_seed",i,"_time.RData",sep="")
  load(file=dir)
  
  l1time_Gaussian <- rbind(l1time_Gaussian,timek_Gaussian)
  l1time_subGaussian <- rbind(l1time_subGaussian,timek_subGaussian)
  
  rm(timek_Gaussian,timek_subGaussian)
}

l1time_Gaussian_mean <- apply(l1time_Gaussian,2,mean)
l1time_subGaussian_mean <- apply(l1time_subGaussian,2,mean)

l1time_Gaussian_sd <- apply(l1time_Gaussian,2,sd)
l1time_subGaussian_sd <- apply(l1time_subGaussian,2,sd)

l1time_Gaussian_mean
l1time_Gaussian_sd

l1time_subGaussian_mean
l1time_subGaussian_sd


