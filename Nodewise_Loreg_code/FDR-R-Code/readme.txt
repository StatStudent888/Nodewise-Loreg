# This document contains 5 files. 
# The adaptive z-value based procedure is implemented in the function "AdaptZ.func". 
# The inputs are the z-values from m tests and a given FDR level. 
# The function will return the estimated local fdr values,  
  # the threshold for the local fdr values 
  # the number of rejected hypotheses
  # and the set of indices of the rejected hypotheses. 
# Here is an example. 


source("epsest.func.R.txt")
source("lin.itp.R.txt")
source("adpt.cutz.R.txt")
source("adaptZ.func.R.txt")

hiv=scan("hivdata.txt")
FDR<-0.05
adaptiveZ<-adaptZ.func(hiv, FDR)
# the threshold
threshold<-adaptiveZ$th
# number of rejected hypotheses
k<-adaptiveZ$nr
# the rejected hypotheses
rh<-adaptiveZ$re