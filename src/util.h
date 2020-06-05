#ifndef UTIL_H
#define UTIL_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' @name util
//' @export
// [[Rcpp::export]]
arma::mat pinv(const arma::mat& x);

//' @name util
//' @export
// [[Rcpp::export]]
double logchoose(const arma::vec& x);

//' @name util
//' @export
// [[Rcpp::export]]
arma::vec normalize(const arma::vec& x, bool na_rm = true);

#endif

