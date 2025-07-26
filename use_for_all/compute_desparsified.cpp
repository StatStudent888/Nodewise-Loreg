// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

// A function to compute M^(1/2) for an input matrix M
// [[Rcpp::export]]
arma::mat matrix_sqrt(const arma::mat& mat) {

  arma::vec eigenvalues;
  arma::mat eigenvectors;

  if (!arma::eig_sym(eigenvalues, eigenvectors, mat)) {
    stop("Eigenvalue decomposition failed");
  }

  arma::vec sqrt_eigenvalues = arma::sqrt(eigenvalues);

  arma::mat sqrt_mat = eigenvectors * arma::diagmat(sqrt_eigenvalues) * eigenvectors.t();
  return sqrt_mat;
}

// A function to compute M^(-1) for an input matrix M
// [[Rcpp::export]]
arma::mat matrix_inverse(const arma::mat& mat) {

  return arma::inv(mat);
}

// A function to compute the true length for confidence interval
// [[Rcpp::export]]
Rcpp::List compute_desparsified(int j, Eigen::MatrixXd omega2, Eigen::MatrixXd Sigma_sqrt, double alpha, int n, int p) {

  Eigen::VectorXd ltrue_Gaussian_desparsified = Eigen::VectorXd::Zero(j + 1);
  Eigen::VectorXd ltrue_subGaussian_desparsified = Eigen::VectorXd::Zero(j + 1);

  double qnorm_const = 2 * R::qnorm(1 - alpha / 2, 0, 1, 1, 0) / std::sqrt(n);

  for (int i = 0; i <= j; ++i) {
    double sigmatrue_Gaussian_ij = std::sqrt(omega2(i, i) * omega2(j, j) + omega2(i, j) * omega2(i, j));
    ltrue_Gaussian_desparsified[i] = qnorm_const * sigmatrue_Gaussian_ij;

    Eigen::MatrixXd temp = 0.5 * (omega2.col(i) * omega2.row(j) + omega2.col(j) * omega2.row(i));
    Eigen::MatrixXd A = Sigma_sqrt * temp * Sigma_sqrt;

    double trace_AA = A.squaredNorm(); 
    double trace_A_diag2 = A.diagonal().squaredNorm(); 

    double sigmatrue_subGaussian_ij = std::sqrt(2 * trace_AA - 1.2 * trace_A_diag2);
    ltrue_subGaussian_desparsified[i] = qnorm_const * sigmatrue_subGaussian_ij;
  }

  return Rcpp::List::create(
      Rcpp::Named("ltrue_Gaussian_desparsified") = ltrue_Gaussian_desparsified,
      Rcpp::Named("ltrue_subGaussian_desparsified") = ltrue_subGaussian_desparsified);
}


