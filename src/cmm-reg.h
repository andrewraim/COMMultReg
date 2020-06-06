#ifndef CMM_REG_H
#define CMM_REG_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' @export
// [[Rcpp::export]]
Rcpp::List gunterize(const arma::umat& X, bool all = false);

//' @export
// [[Rcpp::export]]
Rcpp::List loglik_score_fim_cmm(const Rcpp::List& par,
	const Rcpp::List& dat_xform, unsigned int baseline);

#endif
