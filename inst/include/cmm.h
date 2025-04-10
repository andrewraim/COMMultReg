#ifndef COMMULTREG_CMM_H
#define COMMULTREG_CMM_H

#include "cmb.h"
#include "MultIterator.h"

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace COMMultReg {

/*
* The CMM conditionals are CMB random variables. We can use this fact to
* sample from CMM using a Gibbs sampler.
*/
// [[Rcpp::export]]
inline arma::mat r_cmm_internal(unsigned int n, unsigned int m, const arma::vec& p,
	double nu, unsigned int burn, unsigned int thin, const arma::vec& x_init,
	unsigned int report)
{
	if (!arma::all(p > 0)) {
		Rcpp::stop("All coordinates of p must be positive");
	}

	unsigned int k = p.n_elem;
	const arma::vec& p_norm = p / sum(p);

	unsigned int R = n * thin + burn;
	arma::vec x = x_init;
	arma::mat x_hist(n, k);
	unsigned int idx_keep = 0;

	for (unsigned int r = 0; r < R; r++) {
		if ((r+1) % report == 0) {
			Rprintf("Entering step %d of the MCMC\n", r+1);
		}

		for (unsigned int j = 0; j < k-1; j++) {
			unsigned int ss = 0;
			for (unsigned int l = 0; l < k-1; l++) {
				ss += x(l) * (l != j);
			}
			double m_star = m - ss;
			double p_star = p_norm(j) / (p_norm(j) + p_norm(k-1));
			double nu_star = nu;
			const arma::vec& q = r_cmb(1, m_star, p_star, nu_star);
			x(j) = q(0);
			x(k-1) = m_star - x(j);
		}

		if ((r+1) > burn && ((r+1) % thin == 0)) {
			x_hist.row(idx_keep) = x.t();
			idx_keep++;
		}
	}

	return x_hist;
}

//' @name cmm
//' @export
// [[Rcpp::export]]
inline double d_cmm(const arma::vec& x, const arma::vec& p, double nu,
	bool log = false, bool normalize = true)
{
	size_t k = x.n_elem;
	if (k != p.n_elem) {
		Rcpp::stop("Must have dim(x) == dim(p)");
	}

	size_t m = arma::sum(x);
	double ll = nu*lgamma(m+1) - nu*sum(lgamma(x+1)) + arma::dot(x, arma::log(p));

	// Don't store the entire multinomial sample space; iterate through it.
	// Skip this if we don't need the probabilties to be normalized (e.g. if
	// it can be done more efficiently outside).
	if (normalize) {
		MultIterator itr(k, m);
		double sumprob = 0;
		for (; !itr.is_end(); itr.increment()) {
			// Note: tried this with ivec, but got different (wrong) answers
			const arma::vec& xx = itr.getCounts();
			sumprob += exp(nu*lgamma(m+1) - nu*sum(lgamma(xx+1)) + arma::dot(xx, arma::log(p)));
		}
		ll -= std::log(sumprob);
	}

	return log ? ll : exp(ll);
}

//' @name cmm
//' @export
// [[Rcpp::export]]
inline double normconst_cmm(unsigned int m, const arma::vec& p, double nu,
	bool log = false)
{
	size_t k = p.n_elem;
	arma::vec x = arma::zeros(k);
	x(0) = m;
	double out = d_cmm(x, p, nu, true, false) - d_cmm(x, p, nu, true, true);
	return log ? out : exp(out);
}


//' Density for CMM random sample
//' 
//' Compute individual density contributions for 
//' \deqn{
//' \bm{X}_i \sim \textrm{CMM}_k(m_i, \bm{p}_i, \nu_i),
//' \quad i = 1, \ldots, n.
//' }
//' 
//' @param X An \eqn{n \times k} matrix of outcomes, where the \eqn{i}th row
//' \eqn{\bm{x}_i^\top} represents the \eqn{i}th observation.
//' @param P An \eqn{n \times k} matrix, where the \eqn{i}th row
//' \eqn{\bm{p}_i^\top} represents the probability parameter for the
//' \eqn{i}th observation.
//' @param nu An \eqn{n}-dimensional vector of dispersion parameters
//' \eqn{\nu_1, \ldots, \nu_n}
//' @param log \code{TRUE} or \code{FALSE}; if \code{TRUE}, return the
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
//' according to input arguments. The value of each
//' \eqn{m_i} is assumed to be \eqn{\sum_{j=1}^k x_{ij}}.
//'
//' @details
//' The entire computation for this function is done in C++, and therefore
//' may be more efficient than calling \code{d_cmm} in a loop from R.
//'
//' @examples
//' set.seed(1234)
//' 
//' n = 20
//' m = rep(10, n)
//' k = 3
//' 
//' x = rnorm(n)
//' X = model.matrix(~ x)
//' beta = matrix(NA, 2, k-1)
//' beta[1,] = -1
//' beta[2,] = 1
//' P = t(apply(X %*% beta, 1, inv_mlogit))
//' 
//' w = rnorm(n)
//' W = model.matrix(~ x)
//' gamma = c(1, -0.1)
//' nu = X %*% gamma
//' 
//' y = matrix(NA, n, k)
//' for (i in 1:n) {
//'     y[i,] = r_cmm(1, m[i], P[i,], nu[i], burn = 200)
//' }
//' 
//' d_cmm_sample(y, P, nu, log = TRUE)
//' 
//' @export
// [[Rcpp::export]]
inline arma::vec d_cmm_sample(const arma::mat& X, const arma::mat& P,
	const arma::vec& nu, bool log = false, bool normalize = true)
{
	unsigned int n = X.n_rows;
	unsigned int k = X.n_cols;
	if (n != P.n_rows || k != P.n_cols || n != nu.n_elem) {
		Rcpp::stop("Must have length(nu) == nrow(X) and dim(X) == dim(P)");
	}

	arma::vec out(n);
	for (size_t i = 0; i < n; i++) {
		out(i) = d_cmm(X.row(i).t(), P.row(i).t(), nu(i), log, normalize);
	}

	return out;
}

}

#endif
