adpt.cutz<-function(lfdr, q)
{
  # the input
    # lfdr the vector of local fdr statistics
    # q the desired FDR level
  # the output is a list with
    # the first element (st.lfdr) the sorted local fdr values
    # the second element (k) the number of hypotheses to be rejected
    # the third element (lfdrk) the threshold for the local fdr values
    # the fourth element (reject) the set of indices of the rejected hypotheses
    # the fifth element (accept) the set of indices of the accepted hypotheses
  
  m=length(lfdr)
  st.lfdr<-sort(lfdr)
  k=1
  while(k<m && (1/k)*sum(st.lfdr[1:k])<q){
    k=k+1
  }
  k<-k-1
  lfdrk<-st.lfdr[k]
  reject<-which(lfdr<=lfdrk)
  accept<-which(lfdr>lfdrk)
  y<-list(sf=st.lfdr, nr=k, thr=lfdrk, re=reject, ac=accept)
  return (y)
}
  