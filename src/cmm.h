#ifndef CMM_H
#define CMM_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat r_cmm_internal(unsigned int n, unsigned int m, const arma::vec& p,
	double nu, unsigned int burn, unsigned int thin,
	const arma::vec& x_init, unsigned int report_period);

//' @name cmm
//' @export
// [[Rcpp::export]]
double d_cmm(const arma::vec& x, unsigned int m, const arma::vec& p,
	double nu, bool take_log = false, bool normalize = true);

//' @name cmm
//' @export
// [[Rcpp::export]]
double d_cmm_normconst(unsigned int m, const arma::vec& p, double nu,
	bool take_log = false);

//' @name cmm
//' @export
// [[Rcpp::export]]
arma::vec d_cmm_vectorized(const arma::mat& x, const arma::vec& m, const arma::mat& p,
	const arma::vec& nu, bool take_log = false, bool normalize = true);

//' @name cmm
//' @export
// [[Rcpp::export]]
double loglik_cmm(const arma::mat& y, const arma::vec& m,
	const arma::mat& P_mat, const arma::vec& nu);

#endif

