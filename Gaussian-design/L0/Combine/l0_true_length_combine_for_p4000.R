COL <- 1:4000
n <- 200
p <- 4000

library(BiocParallel)
library(parallel)

ncores=40
mcparam=SnowParam(workers=ncores)

load_function <- function(i,n,p){
  dir <- paste("L0_true_length/L0_Gaussian_subGaussian_band_n",n,"_p",p,"_col",i,"_true_length.RData",sep="")
  load(file=dir)
  
  return(list(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified))
  
  rm(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified)
}

# load_function <- function(i,n,p){
#   dir <- paste("L0_true_length/L0_Gaussian_subGaussian_random_n",n,"_p",p,"_col",i,"_true_length.RData",sep="")
#   load(file=dir)
#   
#   return(list(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified))
#   
#   rm(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified)
# }
# 
# load_function <- function(i,n,p){
#   dir <- paste("L0_true_length/L0_Gaussian_subGaussian_hub_n",n,"_p",p,"_col",i,"_true_length.RData",sep="")
#   load(file=dir)
#   
#   return(list(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified))
#   
#   rm(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified)
# }
# 
# load_function <- function(i,n,p){
#   dir <- paste("L0_true_length/L0_Gaussian_subGaussian_cluster_n",n,"_p",p,"_col",i,"_true_length.RData",sep="")
#   load(file=dir)
#   
#   return(list(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified))
#   
#   rm(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified)
# }

result <- bplapply(COL,load_function,n,p)

ltrue_Gaussian_desparsified_4000 <- numeric(p*(p+1)/2)
ltrue_subGaussian_desparsified_4000 <- numeric(p*(p+1)/2)

current_pos <- 1

for(i in 1:4000){
  ltrue_Gaussian_desparsified <- result[[i]][[1]]
  ltrue_subGaussian_desparsified <- result[[i]][[2]]
  
  len <- length(ltrue_Gaussian_desparsified)
  
  ltrue_Gaussian_desparsified_4000[current_pos:(current_pos + len - 1)] <- ltrue_Gaussian_desparsified
  ltrue_subGaussian_desparsified_4000[current_pos:(current_pos + len - 1)] <- ltrue_subGaussian_desparsified
  
  current_pos <- current_pos + len
  # print(i)
}

ltrue_Gaussian_desparsified <- ltrue_Gaussian_desparsified_4000
ltrue_subGaussian_desparsified <- ltrue_subGaussian_desparsified_4000

Gaussian_matrix <- matrix(0,nrow=p,ncol=p)
Gaussian_matrix[which(upper.tri(Gaussian_matrix, diag = TRUE))] <- ltrue_Gaussian_desparsified
Gaussian_matrix[lower.tri(Gaussian_matrix,diag = FALSE)] <- t(Gaussian_matrix)[lower.tri(Gaussian_matrix,diag = FALSE)]
ltrue_Gaussian_desparsified <- as.vector(Gaussian_matrix)

subGaussian_matrix <- matrix(0,nrow=p,ncol=p)
subGaussian_matrix[which(upper.tri(subGaussian_matrix, diag = TRUE))] <- ltrue_subGaussian_desparsified
subGaussian_matrix[lower.tri(subGaussian_matrix,diag = FALSE)] <- t(subGaussian_matrix)[lower.tri(subGaussian_matrix,diag = FALSE)]
ltrue_subGaussian_desparsified <- as.vector(subGaussian_matrix)

save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_band_n",n,"_p",p,"_normality.RData",sep=""))

# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_random_n",n,"_p",p,"_normality.RData",sep=""))
# 
# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_hub_n",n,"_p",p,"_normality.RData",sep=""))
# 
# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_cluster_n",n,"_p",p,"_normality.RData",sep=""))

n <- 400
ltrue_Gaussian_desparsified <- ltrue_Gaussian_desparsified/sqrt(2)
ltrue_subGaussian_desparsified <- ltrue_subGaussian_desparsified/sqrt(2)

save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_band_n",n,"_p",p,"_normality.RData",sep=""))

# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_random_n",n,"_p",p,"_normality.RData",sep=""))
# 
# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_hub_n",n,"_p",p,"_normality.RData",sep=""))
# 
# save(ltrue_Gaussian_desparsified,ltrue_subGaussian_desparsified,file=paste("L0_desparsified_true_Gaussian_subGaussian_cluster_n",n,"_p",p,"_normality.RData",sep=""))
