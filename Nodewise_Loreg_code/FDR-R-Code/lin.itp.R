lin.itp<-function(x, X, Y){
  ## x: the coordinates of points where the density needs to be interpolated
  ## X: the coordinates of the estimated densities
  ## Y: the values of the estimated densities
  ## the output is the interpolated densities
  x.N<-length(x)
  X.N<-length(X)
  y<-rep(0, x.N)
  for (k in 1:x.N){
    i<-max(which((x[k]-X)>=0))
    if (i<X.N)
      y[k]<-Y[i]+(Y[i+1]-Y[i])/(X[i+1]-X[i])*(x[k]-X[i])
    else 
      y[k]<-Y[i]
  }
  return(y)
}
