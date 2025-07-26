#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// sdar function
// [[Rcpp::export]]
List  sdar(arma::mat X,arma::mat y,int T,int J,int isnorm)
  
{
  
  // int  alpha=0;
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
  // mat e;
  mat G;
  mat ed;
  
  while(isnorm==1)
  {
    mat nX=normalise(X);  //X normalization;n*p matrix 
  }
  mat nX=X;
  
  int n =nX.n_rows;
  int p =nX.n_cols;
  
  
  mat nXT = nX.t();
  mat nXty = nXT * y; 
  
  //initial value
  mat initial= zeros(p,1);
  S=linspace<uvec>(0, p-1, p);
  mat ebeta= zeros(p,1);
  // mat pd=initial+ nXty-nXT*(nX*initial);
  mat pd=nXty;
  // mat a = sort(abs(pd),"descend");
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
    
    S=linspace<uvec>(0, p-1, p);
    for(int i=0;i<T;i++)//XtyAc = Xty(Ac)
    {
      int c=Ac[i];
      nXtyAc.row(i)=nXty.row(c);
      S.row(c)=p+10;
      nXac.col(i)=nX.col(c);
      
    }//XtyAc = Xty(Ac)
    Ic=find(S<p+10);
    
    
    nXacT = nXac.t();
    G=nXacT*nXac;
    
    tbetaac=solve(G,nXtyAc);
    td=nXty-nXT*(nXac*tbetaac);
    
    //ebeta
    for(int i=0;i<T;i++)//ebeta(Ac) = tbetaac;
    {
      ebeta(Ac(i)) = tbetaac(i);
    }
    
    
    for(int i=0;i<p-T;i++) //tdic = td(Ic);
      {
        ed(Ic(i)) = td(Ic(i));
      }
    
    pd =ed +ebeta;
    
    b = sort_index(abs(pd),"descend");
    oAc=Ac;
    Ac = b.rows(0,T-1);
    stop= approx_equal(Ac,oAc, "absdiff", 0.0);
  }
  
  
  return List::create(n,p,ebeta,ed,nlter,stop,Ac);
  
}

// asdar function
// [[Rcpp::export]]
arma::vec asdar(const arma::mat& X, const arma::mat& nX,
                const arma::vec& dx, const arma::vec& y,
                int kmax, int n, int p) {
    arma::mat Newbeta(p - 1, kmax + 1, arma::fill::zeros); 
    arma::vec HBic(kmax + 1, arma::fill::zeros); 

    for (int k = 0; k <= kmax; ++k) {
        if (k == 0) {
            arma::vec ebeta = arma::zeros(p - 1);
            Newbeta.col(k) = ebeta;

            arma::vec residual = y;
            double residual_norm = std::pow(arma::norm(residual, 2), 2) / n;
            double hbick = std::log(residual_norm) + std::log(std::log(n)) * std::log(p - 1) * k / n;
            HBic[k] = hbick;
        } else {
            List a = sdar(nX, y, k, 20, 0);
            arma::vec a_vector = Rcpp::as<arma::vec>(a[2]);
            arma::vec ebeta = dx % a_vector;

            Newbeta.col(k) = ebeta;

            arma::vec residual = X * ebeta - y;
            double residual_norm = std::pow(arma::norm(residual, 2), 2) / n;
            double hbick = std::log(residual_norm) + std::log(std::log(n)) * std::log(p - 1) * k / n;
            HBic[k] = hbick;
        }
    }

    int min_index = std::distance(HBic.begin(), std::min_element(HBic.begin(), HBic.end()));
    return Newbeta.col(min_index);
}

// asdar_iter function
// [[Rcpp::export]]
List asdar_iter(const arma::mat& X, const arma::mat& nX,
                const arma::vec& dx, const arma::vec& y,
                int kmax, int n, int p) {
    arma::mat Newbeta(p - 1, kmax + 1, arma::fill::zeros); 
    arma::vec HBic(kmax + 1, arma::fill::zeros); 
    arma::ivec iter_count(kmax + 1, arma::fill::zeros);
    arma::ivec stop_flags(kmax + 1, arma::fill::zeros);

    for (int k = 0; k <= kmax; ++k) {
        if (k == 0) {
            arma::vec ebeta = arma::zeros(p - 1);
            Newbeta.col(k) = ebeta;

            arma::vec residual = y;
            double residual_norm = std::pow(arma::norm(residual, 2), 2) / n;
            double hbick = std::log(residual_norm) + std::log(std::log(n)) * std::log(p - 1) * k / n;
            HBic[k] = hbick;
            iter_count[k] = 0;
            stop_flags[k] = 1;
        } else {
            List a = sdar(nX, y, k, 20, 0);
            arma::vec a_vector = Rcpp::as<arma::vec>(a[2]);
            int niter = Rcpp::as<int>(a[4]);
            bool stop_flag = Rcpp::as<bool>(a[5]);
            arma::vec ebeta = dx % a_vector;

            Newbeta.col(k) = ebeta;

            arma::vec residual = X * ebeta - y;
            double residual_norm = std::pow(arma::norm(residual, 2), 2) / n;
            double hbick = std::log(residual_norm) + std::log(std::log(n)) * std::log(p - 1) * k / n;
            HBic[k] = hbick;
            iter_count[k] = niter;
            stop_flags[k] = stop_flag ? 1 : 0;
        }
    }

    int min_index = std::distance(HBic.begin(), std::min_element(HBic.begin(), HBic.end()));
    int iter_flag;
    if (iter_count[min_index] < 20) {
        iter_flag = 1;  
    } else if (iter_count[min_index] == 20 && stop_flags[min_index] == 1) {
        iter_flag = 1;  
    } else {
        iter_flag = 0;  
    }

    return Rcpp::List::create(
        Rcpp::Named("beta") = Newbeta.col(min_index),
        Rcpp::Named("iter_count") = iter_count[min_index],
        Rcpp::Named("iter_flag") = iter_flag
    );
}

// Estimate precision matrix
// [[Rcpp::export]]
List process(const arma::mat& data,
                      const arma::mat& data_normalize,
                      const arma::vec& diag,
                      int n, int p) {
    arma::vec sigma(p, arma::fill::zeros); 
    arma::mat sigma_c(p - 1, p, arma::fill::zeros); 

    for (int j = 1; j <= p; ++j) {
        arma::mat X = data;
        X.shed_col(j - 1);
        arma::vec y = data.col(j - 1);
        arma::mat nX = data_normalize;
        nX.shed_col(j - 1);
        arma::vec dx = diag;
        dx.shed_row(j - 1);

        arma::vec best_ebeta = asdar(X, nX, dx, y, 20, n, p);
        double sigmak = n / std::pow(arma::norm(y - X * best_ebeta, 2), 2);
        arma::vec sigma_ck = -sigmak * best_ebeta;
        
        sigma[j - 1] = sigmak; 
        sigma_c.col(j - 1) = sigma_ck; 
    }

    return List::create(Named("sigma") = sigma,
                        Named("sigma_c") = sigma_c);
}

// A function to compute M^(-1) for an input matrix M
// [[Rcpp::export]]
arma::mat matrix_inverse(const arma::mat& M) {
    arma::mat result = arma::inv(M);
    return result;
}

// A function to compute M^(1/2) for an input matrix M
// [[Rcpp::export]]
arma::mat matrix_sqrt(const arma::mat& M) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, M);
    arma::mat D = arma::diagmat(arma::sqrt(eigval));
    return eigvec * D * eigvec.t();
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec computeZValueLowGaussian(const arma::mat& sdar_T, const arma::mat& sdar_F_sigma3, int n, int p) {
    int num_elements = (p * (p - 1)) / 2;  
    
    arma::vec zvalue_low(num_elements);
    int index = 0;  

    for (int j = 0; j < (p - 1); ++j) {
        for (int i = (j + 1); i < p; ++i) {
            double T_hat_ij = sdar_T(i, j);

            double sigmahat_ij2 = sdar_F_sigma3(i, i) * sdar_F_sigma3(j, j) + 
                                  std::pow(sdar_F_sigma3(i, j), 2);

            double sigmahat_ji2 = sdar_F_sigma3(i, i) * sdar_F_sigma3(j, j) + 
                                  std::pow(sdar_F_sigma3(j, i), 2);

            double sigmahat_ij_lowertri = std::sqrt((sigmahat_ij2 + sigmahat_ji2) / 2);

            double zvalue_low_ij = std::sqrt(n) * T_hat_ij / sigmahat_ij_lowertri;

            zvalue_low(index++) = zvalue_low_ij;
        }
    }

    return zvalue_low;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat matrixMultiply(const arma::mat& A, const arma::mat& B) {
    arma::mat result = A * B;

    return result;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec computeZValueDesparsifiedGaussian(const arma::mat& sdar_F_lowertri1, 
                             const arma::mat& sdar_T, 
                             const arma::mat& sdar_F_sigma3, 
                             const arma::vec& sigma, 
                             int n, 
                             int p) {
  std::vector<double> zvalue_desparsified_min;

  for (int j = 0; j < p; ++j) {
    // Exclude row j and extract column j from sdar_F_lowertri1
    arma::vec sdar_F_lowertri1_col_j = sdar_F_lowertri1.col(j);
    sdar_F_lowertri1_col_j.shed_row(j);

    // Get the non-zero indices (Ajhat)
    arma::uvec Ajhat = find(sdar_F_lowertri1_col_j != 0);
    int nj = Ajhat.n_elem;
    if (nj == 0) continue;

    for (int i = 0; i < nj; ++i) {
      int idx = Ajhat[i];

      // Extract necessary elements
      arma::vec sdar_T_col_j = sdar_T.col(j);
      sdar_T_col_j.shed_row(j);
      double T_hat_ij = sdar_T_col_j(idx);

      arma::vec sdar_F_sigma3_col_j = sdar_F_sigma3.col(j);
      sdar_F_sigma3_col_j.shed_row(j);
      double theta_hat_ij = sdar_F_sigma3_col_j(idx);

      arma::vec sdar_F_sigma3_row_j = sdar_F_sigma3.row(j).t();
      sdar_F_sigma3_row_j.shed_row(j);
      double theta_hat_ji = sdar_F_sigma3_row_j(idx);

      // Exclude row j and column j from sdar_F_sigma3
      arma::mat sdar_F_sigma3_excl = sdar_F_sigma3;
      sdar_F_sigma3_excl.shed_row(j);
      sdar_F_sigma3_excl.shed_col(j);

      // Calculate sigmahat values
      double sigmahat_ij2 = sdar_F_sigma3_excl(idx, idx) * sigma[j] + std::pow(theta_hat_ij, 2);
      double sigmahat_ji2 = sdar_F_sigma3_excl(idx, idx) * sigma[j] + std::pow(theta_hat_ji, 2);
      double sigmahat_ij_lowertri = std::sqrt((sigmahat_ij2 + sigmahat_ji2) / 2);

      // Calculate z-value
      double zvalue_ij = std::sqrt(n) * T_hat_ij / sigmahat_ij_lowertri;

      // Append to result vector
      zvalue_desparsified_min.push_back(zvalue_ij);
    }
  }

  return arma::vec(zvalue_desparsified_min);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec computeZValueLowSubGaussian(const arma::mat& sdar_T, 
                             const arma::mat& sdar_F_sigma3, 
                             const arma::mat& dataSigma2, 
                             int n, int p) {
  
  int num_elements = p * (p - 1) / 2; 
  arma::vec zvalue_low(num_elements);

  int index = 0; 
  
  for (int j = 0; j < (p - 1); ++j) {
    for (int i = (j + 1); i < p; ++i) {
      double T_hat_ij = sdar_T(i, j);

      double term1 = arma::dot(dataSigma2.col(i),dataSigma2.col(j)) / n;
      double term2 = 0.5 * (std::pow(sdar_F_sigma3(i, j), 2) + std::pow(sdar_F_sigma3(j, i), 2));
      double sigmahat_ij2 = term1 - term2;

      double sigmahat_ij_lowertri = std::sqrt(sigmahat_ij2);

      double zvalue_low_ij = std::sqrt(n) * T_hat_ij / sigmahat_ij_lowertri;

      zvalue_low(index++) = zvalue_low_ij;
    }
  }

  return zvalue_low;
}


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List desparsified_confidence_intervals_Gaussian(const arma::mat& sdar_T, 
                                             const arma::mat& sdar_F_sigma, 
                                             const arma::mat& omega2, 
                                             const int n, 
                                             const int p, 
                                             const double alpha) {
  
    int num_elements = (p * (p + 1)) / 2; 
    
    arma::vec z_desparsified(num_elements);
    arma::vec l_desparsified(num_elements);
    arma::vec cov_desparsified(num_elements);

    int zindex = 0;
    int lindex = 0;
    int covindex = 0;
    
    for (int j = 0; j < p; ++j) { 
        for (int i = 0; i <= j; ++i) { 
            
            double sigmahat_ij2 = sdar_F_sigma(i, i) * sdar_F_sigma(j, j) + 
                           0.5 * (std::pow(sdar_F_sigma(i, j), 2) + std::pow(sdar_F_sigma(j, i), 2));
            double sigmahat_ij = std::sqrt(sigmahat_ij2);
            
            double z_ij = std::sqrt(n) * (sdar_T(i, j) - omega2(i, j)) / sigmahat_ij;
            z_desparsified(zindex++) = z_ij;
            
            double qnorm = R::qnorm(1 - (alpha / 2), 0.0, 1.0, true, false);
            double l_ij = 2 * qnorm * sigmahat_ij / std::sqrt(n);
            l_desparsified(lindex++) = l_ij;
            
            double CI_ij_low = sdar_T(i, j) - qnorm * sigmahat_ij / std::sqrt(n);
            double CI_ij_up = sdar_T(i, j) + qnorm * sigmahat_ij / std::sqrt(n);
            
            int cov_ij = (omega2(i, j) >= CI_ij_low && omega2(i, j) <= CI_ij_up) ? 1 : 0;
            cov_desparsified(covindex++) = cov_ij;
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("z_desparsified") = z_desparsified,
                              Rcpp::Named("l_desparsified") = l_desparsified,
                              Rcpp::Named("cov_desparsified") = cov_desparsified);
}


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List desparsified_confidence_intervals_SubGaussian(const arma::mat& sdar_T, 
                                             const arma::mat& sdar_F_sigma, 
                                             const arma::mat& omega2, 
                                             const arma::mat& dataSigma2,
                                             const int n, 
                                             const int p, 
                                             const double alpha) {
  
    int num_elements = (p * (p + 1)) / 2; 
    
    arma::vec z_desparsified(num_elements);
    arma::vec l_desparsified(num_elements);
    arma::vec cov_desparsified(num_elements);

    int zindex = 0;
    int lindex = 0;
    int covindex = 0;
    
    for (int j = 0; j < p; ++j) { 
        for (int i = 0; i <= j; ++i) { 
            double T_hat_ij = sdar_T(i, j);

            double term1 = arma::dot(dataSigma2.col(i),dataSigma2.col(j)) / n;
            double term2 = 0.5 * (std::pow(sdar_F_sigma(i, j), 2) + std::pow(sdar_F_sigma(j, i), 2));
            double sigmahat_ij2 = term1 - term2;
            double sigmahat_ij = std::sqrt(sigmahat_ij2);
            
            double z_ij = std::sqrt(n) * (T_hat_ij - omega2(i, j)) / sigmahat_ij;
            z_desparsified(zindex++) = z_ij;
            
            double qnorm = R::qnorm(1 - (alpha / 2), 0.0, 1.0, true, false);
            double l_ij = 2 * qnorm * sigmahat_ij / std::sqrt(n);
            l_desparsified(lindex++) = l_ij;
            
            double CI_ij_low = T_hat_ij - qnorm * sigmahat_ij / std::sqrt(n);
            double CI_ij_up = T_hat_ij + qnorm * sigmahat_ij / std::sqrt(n);
            
            int cov_ij = (omega2(i, j) >= CI_ij_low && omega2(i, j) <= CI_ij_up) ? 1 : 0;
            cov_desparsified(covindex++) = cov_ij;
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("z_desparsified") = z_desparsified,
                              Rcpp::Named("l_desparsified") = l_desparsified,
                              Rcpp::Named("cov_desparsified") = cov_desparsified);
}


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List L1_desparsified_confidence_intervals_SubGaussian(const arma::mat& NS_T, const arma::mat& NS_F_sigma, const arma::mat& omega2,
 const arma::mat& dataSigma2, int n, int p, double alpha) {

  int num_elements = (p * (p + 1)) / 2;  // Total elements for upper triangular matrix
  arma::vec z_desparsified(num_elements);
  arma::vec l_desparsified(num_elements);
  arma::vec cov_desparsified(num_elements);

  int zindex = 0, lindex = 0, covindex = 0;

  for (int j = 0; j < p; ++j) {
    for (int i = 0; i <= j; ++i) {
      double T_hat_ij = NS_T(i, j);
      double right = arma::dot(dataSigma2.col(i), dataSigma2.col(j))/ n;
      double sigmahat_ij2 = right
                            - 0.5 * (std::pow(NS_F_sigma(i, j), 2) + std::pow(NS_F_sigma(j, i), 2));
      if (sigmahat_ij2 <= 0) {
        double sigmahat_ij2_fix = right
                                  - std::min(std::pow(NS_F_sigma(i, j), 2), std::pow(NS_F_sigma(j, i), 2));
        if (sigmahat_ij2_fix <= 0) {
          z_desparsified(zindex++) = NA_REAL;
          l_desparsified(lindex++) = NA_REAL;
          cov_desparsified(covindex++) = NA_REAL;
        } else {
          double sigmahat_ij = std::sqrt(sigmahat_ij2_fix);
          double z_ij = std::sqrt(n) * (T_hat_ij - omega2(i, j)) / sigmahat_ij;
          z_desparsified(zindex++) = z_ij;
          double l_ij = 2 * R::qnorm(1 - alpha / 2, 0, 1, 1, 0) * sigmahat_ij / std::sqrt(n);
          l_desparsified(lindex++) = l_ij;
          double CI_low = T_hat_ij - R::qnorm(1 - alpha / 2, 0, 1, 1, 0) * sigmahat_ij / std::sqrt(n);
          double CI_up = T_hat_ij + R::qnorm(1 - alpha / 2, 0, 1, 1, 0) * sigmahat_ij / std::sqrt(n);
          cov_desparsified(covindex++) = (omega2(i, j) <= CI_up && omega2(i, j) >= CI_low) ? 1.0 : 0.0;
        }
      } else {
        double sigmahat_ij = std::sqrt(sigmahat_ij2);
        double z_ij = std::sqrt(n) * (T_hat_ij - omega2(i, j)) / sigmahat_ij;
        z_desparsified(zindex++) = z_ij;
        double l_ij = 2 * R::qnorm(1 - alpha / 2, 0, 1, 1, 0) * sigmahat_ij / std::sqrt(n);
        l_desparsified(lindex++) = l_ij;
        double CI_low = T_hat_ij - R::qnorm(1 - alpha / 2, 0, 1, 1, 0) * sigmahat_ij / std::sqrt(n);
        double CI_up = T_hat_ij + R::qnorm(1 - alpha / 2, 0, 1, 1, 0) * sigmahat_ij / std::sqrt(n);
        cov_desparsified(covindex++) = (omega2(i, j) <= CI_up && omega2(i, j) >= CI_low) ? 1.0 : 0.0;
      }
    }
  }

  return List::create(Named("z_desparsified") = z_desparsified,
                      Named("cov_desparsified") = cov_desparsified,
                      Named("l_desparsified") = l_desparsified);
}