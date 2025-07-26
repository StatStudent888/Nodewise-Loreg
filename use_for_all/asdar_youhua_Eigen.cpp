// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <queue>
#include <utility>

using namespace Rcpp;
using namespace arma;

// A function that uses a min-heap sorting algorithm to return the indices of the top T largest values in an input vector
Eigen::VectorXi topTIndices(const Eigen::VectorXd& vec, size_t T) {
    // Priority queue to store pairs of (value, index) sorted by value
    using Pair = std::pair<double, size_t>;
    auto compare = [](const Pair& a, const Pair& b) { return a.first > b.first; };
    std::priority_queue<Pair, std::vector<Pair>, decltype(compare)> minHeap(compare);

    // Iterate through the vector
    for (size_t i = 0; i < vec.size(); ++i) {
        if (minHeap.size() < T) {
            minHeap.emplace(vec[i], i);
        } else if (vec[i] > minHeap.top().first) {
            minHeap.pop();
            minHeap.emplace(vec[i], i);
        }
    }

    // Extract indices from the priority queue
    Eigen::VectorXi indices(T);
    for (size_t i = 0; i < T; ++i) {
        indices(T - i - 1) = minHeap.top().second; // Store in reverse order
        minHeap.pop();
    }

    return indices;
}

// sdar funtion in Eigen version, which is faster than armadillo
Eigen::VectorXd  sdar(const Eigen::MatrixXd& X, const Eigen::MatrixXd& y,int T,int J,int isnorm)
  
{
  
  int  nlter=0;
  bool stop=0;
  
  Eigen::VectorXi Ac, oAc, Ic, S;
  Eigen::VectorXd ed, ebeta, tbetaac, td;
  Eigen::MatrixXd nXacT;
  Eigen::MatrixXd G;
  
  Eigen::MatrixXd nX = (isnorm == 1) ? X.colwise().normalized() : X;
  
  int n =nX.rows();
  int p =nX.cols();
  
  Eigen::MatrixXd nXT = nX.transpose();
  Eigen::VectorXd nXty = nXT * y;
  
  size_t T_size = T;
  //initial value
  Eigen::VectorXd pd = nXty;
  // b = sort_index(abs(pd),"descend");
  Ac = topTIndices(pd.cwiseAbs(), T_size);
  Eigen::VectorXd nXtyAc = Eigen::VectorXd::Zero(T);
  Eigen::MatrixXd nXac = Eigen::MatrixXd::Zero(n, T);
  
  
 
  //loop
  while(stop==0&&nlter<J)
  {
    nlter++;
    ed = Eigen::VectorXd::Zero(p);
    ebeta = Eigen::VectorXd::Zero(p);
    
    S = Eigen::VectorXi::LinSpaced(p, 0, p - 1);
    nXtyAc = nXty(Ac);
    for (int i = 0; i < Ac.size(); ++i) {
            S(Ac[i]) = p + 10;
        }
    nXac = nX(Eigen::all, Ac);
    Ic.resize(p - T); 
    int count = 0;

    for (int i = 0; i < S.size(); ++i) {
        if (S[i] < p + 10) {
            Ic[count] = i; 
            ++count;       
        }
    }
    
    G=nXac.transpose() * nXac;
    
    tbetaac=G.ldlt().solve(nXtyAc);
    td=nXty-nXT*(nXac*tbetaac);

    for (int i = 0; i < Ac.size(); ++i) {
        ebeta(Ac[i]) = tbetaac[i];
    }
    for (int i = 0; i < Ic.size(); ++i) {
        ed(Ic[i]) = td(Ic[i]);
    }
    
    pd =ed +ebeta;
    
    oAc = Ac;
    Ac = topTIndices(pd.cwiseAbs(), T_size);
    stop= ((Ac - oAc).cwiseAbs().sum() == 0);
  }
  
  
  return ebeta;
  
}

// asdar function in Eigen version
// [[Rcpp::export]]
Eigen::VectorXd asdar(const Eigen::MatrixXd& nX, const Eigen::VectorXd& y,
                int kmax, int n, int p) {
    Eigen::MatrixXd Newbeta = Eigen::MatrixXd::Zero(p - 1, kmax + 1);
    Eigen::VectorXd HBic = Eigen::VectorXd::Zero(kmax + 1);
    
    double log_log_n_p = std::log(std::log(n)) * std::log(p - 1) / n;

    for (int k = 0; k <= kmax; ++k) {
        if (k == 0) {
            Newbeta.col(k).setZero();

            double residual_norm = std::pow(y.norm(), 2) / n;
            double hbick = std::log(residual_norm) + log_log_n_p * k ;
            HBic[k] = hbick;
        } else {
            // List a = sdar(nX, y, k, 20, 0);
            Eigen::VectorXd ebeta = sdar(nX, y, k, 20, 0);
            // arma::vec ebeta = dx % a_vector;

            Newbeta.col(k) = ebeta;

            Eigen::VectorXd residual = nX * ebeta - y;
            double residual_norm = std::pow(residual.norm(), 2) / n;
            double hbick = std::log(residual_norm) + log_log_n_p * k;
            HBic[k] = hbick;
        }
    }
    
    int min_index;
    HBic.minCoeff(&min_index);

    // arma::vec mybeta = dx % Newbeta.col(min_index);
    return Newbeta.col(min_index);
}

// Estimate the precision matrix using Eigen
// [[Rcpp::export]]
Rcpp::List process(const Eigen::MatrixXd& data,
                      const Eigen::MatrixXd& data_normalize,
                      const Eigen::VectorXd& diag,
                      int n, int p) {
    Eigen::VectorXd sigma(p);
    Eigen::MatrixXd sigma_c(p - 1, p);

    Eigen::MatrixXd nX(data_normalize.rows(), data_normalize.cols() - 1);
    Eigen::VectorXd dx(diag.size() - 1);

    for (int j = 0; j < p; ++j) {
        Eigen::VectorXd y = data.col(j);

        nX << data_normalize.leftCols(j), data_normalize.rightCols(data_normalize.cols() - j - 1);

        dx << diag.head(j), diag.tail(diag.size() - j - 1);

        Eigen::VectorXd best_ebeta = asdar(nX, y, 20, n, p);
        double sigmak = n / std::pow((y - nX * best_ebeta).norm(), 2);
        Eigen::VectorXd sigma_ck = -sigmak * dx.cwiseProduct(best_ebeta);
        
        sigma(j) = sigmak; 
        sigma_c.col(j) = sigma_ck; 
    }

    return Rcpp::List::create(Rcpp::Named("sigma") = sigma,
                        Rcpp::Named("sigma_c") = sigma_c);
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

