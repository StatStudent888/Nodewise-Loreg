#Generate band graph of precision matrix
####################################
#p-----------Dimension
#rou0,rou1,rou2-----Some constant
##################################
#Generate band matrix
omega_generator1 <- function(p,rou0,rou1,rou2){
  omega <- matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){
        omega[i,j] <- rou0
      }
      else if(abs(i-j)==1){
        omega[i,j] <- rou1
      }
      else if(abs(i-j)==2){
        omega[i,j] <- rou2
      }
      else{omega[i,j] <- 0}
    }
  }
  return(omega)
}