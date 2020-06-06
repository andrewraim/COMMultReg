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
double normconst_cmb(unsigned int m, double p, double nu,
	bool take_log = false);

//' Density for CMB random sample
//' 
//' Compute individual density contributions for
//' \deqn{
//' X_i \sim \textrm{CMB}(m_i, p_i, \nu_i),
//' \quad i = 1, \ldots, n.
//' }
//' 
//' @param x An \eqn{n}-dimensional vector of outcomes
//' @param m An \eqn{n}-dimensional vector \eqn{m_1, \ldots, m_n}
//' @param p An \eqn{n}-dimensional vector of probability parameters
//' \eqn{p_1, \ldots, p_n}
//' @param nu An \eqn{n}-dimensional vector of dispersion parameters
//' \eqn{\nu_1, \ldots, \nu_n}
//' @param take_log \code{TRUE} or \code{FALSE}; if \code{TRUE}, return the
//' value on the log-scale.
//' @param normalize \code{TRUE} or \code{FALSE}; if \code{FALSE}, do not
//' compute or apply the normalizing constant to each density value.
//' 
//' @return
//' A vector of density values
//' \eqn{f(x_1 \mid m_1, p_1, \nu_1), \ldots, f(x_n \mid m_n, p_n, \nu_n),}
//' which may be on the log-scale and/or unnormalized
//' according to input arguments. See \link{cmb}.
//'
//' @examples
//' stop("TBD")
//' 
//' @export
// [[Rcpp::export]]
arma::vec d_cmb_sample(const arma::vec& x, const arma::vec& m,
	const arma::vec& p, const arma::vec& nu);

#endif
