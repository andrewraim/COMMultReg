#ifndef CMB_H
#define CMB_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' @name cmb
//' @export
// [[Rcpp::export]]
arma::vec r_cmb(unsigned int n, unsigned int m, double p, double nu);

//' @name cmb
//' @export
// [[Rcpp::export]]
double d_cmb(unsigned int x, unsigned int m, double p,
	double nu, bool take_log = false, bool normalize = true);

//' @name cmb
//' @export
// [[Rcpp::export]]
double d_cmb_normconst(unsigned int m, double p, double nu,
	bool take_log = false);

//' @name cmb
//' @export
// [[Rcpp::export]]
double loglik_cmb(const arma::vec& y, const arma::vec& m,
	const arma::vec& p, const arma::vec& nu);

#endif
