seed <- scan("seed.txt",what=integer())
n <- 200 #400
p <- 200 #400,1000,4000

l0itercount_median <- {}
l0itercount_IQR <- {}
l0iterflag <- {}

for(i in seed){
  dir <- paste("L0_Gaussian_NumOfIter/L0_Gaussian_band_n",n,"_p",p,"_seed",i,"_NumOfIter.RData",sep="")
  load(file=dir)
  
  count_IQR <- IQR(iter_count)
  l0itercount_median <- rbind(l0itercount_median,count_median)
  l0itercount_IQR <- rbind(l0itercount_IQR,count_IQR)
  l0iterflag <- rbind(l0iterflag,flag_mean)
  
  rm(count_mean,count_sd,count_median,flag_mean,iter_count)
}

l0itercount_median_mean <- mean(l0itercount_median)
l0itercount_IQR_mean <- mean(l0itercount_IQR)
l0iterflag_mean <- mean(l0iterflag)

l0itercount_median_sd <- sd(l0itercount_median)
l0itercount_IQR_sd <- sd(l0itercount_IQR)
l0iterflag_sd <- sd(l0iterflag)

x1 <- paste(sprintf("%0.3f",l0itercount_median_mean),"(",sprintf("%0.3f",l0itercount_median_sd),")",sep="")
x2 <- paste(sprintf("%0.3f",l0itercount_IQR_mean),"(",sprintf("%0.3f",l0itercount_IQR_sd),")",sep="")
x3 <- paste(sprintf("%0.3f",l0iterflag_mean),"(",sprintf("%0.3f",l0iterflag_sd),")",sep="")

x1
x2
x3


