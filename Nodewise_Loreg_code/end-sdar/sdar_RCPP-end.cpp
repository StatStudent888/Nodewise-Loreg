//[[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]

List  sdar(arma::mat X,arma::mat y,int T,int J,int isnorm)
  
{
  
  int  alpha=0;
  int  nlter=0;
  bool stop=0;
  
  uvec Ic;
  uvec oAc;
  uvec Ac;
  uvec b;
  uvec S;
  mat td;
  mat tbetaac;
  mat nXacT;
  mat e;
  mat G;
  mat ed;
  
  while(isnorm==1)
  {
    mat nX=normalise(X);  //X normalization;n*p matrix 
  }
  mat nX=X;
  
  int n =nX.n_rows;
  int p =nX.n_cols;
  //cout<<"n="<<n<<endl;
  //cout<<"p="<<p<<endl;
  
  
  mat nXT = nX.t();
  mat nXty = nXT * y; 
  
  //initial value
  mat initial= zeros(p,1);
  S=linspace<uvec>(0, p-1, p);
  mat ebeta= zeros(p,1);
  mat pd=initial+ nXty-nXT*(nX*initial);
  mat a = sort(abs(pd),"descend");
  b = sort_index(abs(pd),"descend");
  Ac = b.rows(0,T-1);
  mat nXtyAc=zeros(T,1);
  mat tdic=zeros(p-T,1);
  mat nXac=zeros(n,T);
  
  
 
  //loop
  while(stop==0&&nlter<J)
  {
    nlter++;
    ed= zeros(p,1);
    ebeta= zeros(p,1);
    
    for(int i=0;i<T;i++)//��XtyAc = Xty(Ac)��
    {
      int c=Ac[i];
      nXtyAc.row(i)=nXty.row(c);
      
    }//��XtyAc = Xty(Ac)��
    
    S=linspace<uvec>(0, p-1, p);
    //S.print();
    for(int i=0;i<T;i++)
    {
      int c=Ac[i];
      S.row(c)=p+10;
    }
    Ic=find(S<p+10);
    
    
    for(int i=0;i<T;i++)//��Xac = X(:,Ac);��
    {
      int c=Ac[i];
      nXac.col(i)=nX.col(c);
      
    }//��Xac = X(:,Ac);��
    
    nXacT = nXac.t();
    e=eye<mat>(T,T);
    G=nXacT*nXac+alpha*e;
    
    tbetaac=solve(G,nXtyAc);
    
    //ebeta
    for(int i=0;i<T;i++)//��ebeta(Ac) = tbetaac;��
    {
      int c=Ac[i];//
      ebeta.row(c)=tbetaac.row(i);
    }//
    
    td=nXty-nXT*(nXac*tbetaac);
    
    for(int i=0;i<p-T;i++) //tdic = td(Ic);
      {
      int c=Ic[i];
      tdic.row(i)=td.row(c);
      }//tdic = td(Ic);
    
    for(int i=0;i<p-T;i++)//ed(Ic) = tdic;
      {
      int c=Ic[i];
      ed.row(c)=tdic.row(i);
      }//ed(Ic) = tdic;
    
    pd =ed +ebeta;
    
    b = sort_index(abs(pd),"descend");
    oAc=Ac;
    Ac = b.rows(0,T-1);
    oAc=sort(abs(oAc));
    Ac=sort(abs(Ac));
    stop= approx_equal(Ac,oAc, "absdiff", 0.0);
    
    //oAc.print("oAC=");
    //Ac.print("AC=");
    //cout << " stop= " <<stop<< endl;
    //cout<<"nlter="<<nlter<<endl;
  }
  
  
  return List::create(n,p,ebeta,ed,nlter,Ac);
  
}//main function
