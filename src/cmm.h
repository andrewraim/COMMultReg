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
double normconst_cmm(unsigned int m, const arma::vec& p, double nu,
	bool take_log = false);

//' Density for CMM random sample
//' 
//' Compute individual density contributions for observations in an independent, but not necessarily
//' identically distributed, CMM random sample.
//' 
//' @param X An \eqn{n \times k} matrix of outcomes, where the \eqn{i}th row
//' \eqn{\bm{x}_i^\top} represents the \eqn{i}th observation.
//' @param m An \eqn{n}-dimensional vector \eqn{m_1, \ldots, m_n}
//' @param P An \eqn{n \times k} matrix, where the \eqn{i}th row
//' \eqn{\bm{p}_i^\top} represents the probability parameter for the
//' \eqn{i}th observation.
//' @param nu An \eqn{n}-dimensional vector of dispersion parameters
//' \eqn{\nu_1, \ldots, \nu_n}
//' @param take_log \code{TRUE} or \code{FALSE}; if \code{TRUE}, return the
//' value on the log-scale.
//' @param normalize \code{TRUE} or \code{FALSE}; if \code{FALSE}, do not
//' compute or apply the normalizing constant to each density value.
//' 
//' @return
//' A vector of density values
//' \eqn{
//' f(\bm{x}_1^\top \mid m_1, \bm{p}_1^\top, \nu_1),
//' \ldots,
//' f(\bm{x}_n^\top \mid m_n, \bm{p}_n^\top, \nu_n),
//' }
//' which may be on the log-scale and/or unnormalized
//' according to input arguments.
//'
//' @examples
//' stop("TBD")
//' 
// [[Rcpp::export]]
arma::vec d_cmm_sample(const arma::mat& X, const arma::vec& m, const arma::mat& P,
	const arma::vec& nu, bool take_log = false, bool normalize = true);

#endif

