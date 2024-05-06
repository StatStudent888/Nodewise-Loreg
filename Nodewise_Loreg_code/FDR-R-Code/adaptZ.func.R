adaptZ.func<-function(zv, q)
{
  # the input
    # zv is the z-values transformed from m tests
    # q is the desired FDR level
  # the output is a list with
    # the first element (st.lfdr) the sorted local fdr values
    # the second element (k) the number of hypotheses to be rejected
    # the third element (lfdrk) the threshold for the local fdr values
    # the fourth element (reject) the set of indices of the rejected hypotheses
    # the fifth element (accept) the set of indices of the accepted hypotheses
   ## the estimates for the local fdr statistics
  # density estimates 
  zv.ds<-density(zv, from=min(zv)-10, to=max(zv)+10, n=1000)
  # linear interpolation
  zv.ds<-lin.itp(zv, zv.ds$x, zv.ds$y)
  mu<-0
  s<-1
  zv.p0<-1-epsest.func(zv, mu, s)
  zv.lfdr<-zv.p0*dnorm(zv, mu, s)/zv.ds
  y<-adpt.cutz(zv.lfdr, q)
  return (list(threshold=y$thr,lfdr=zv.lfdr,ac=y$ac))
}
  
