#ifndef UTIL_H
#define UTIL_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat pinv(const arma::mat& x);

// [[Rcpp::export]]
double logchoose(const arma::vec& x);

// [[Rcpp::export]]
arma::vec normalize(const arma::vec& x, bool na_rm = true);

#endif

