#ifndef CMM_H
#define CMM_H

/*
* The CMM conditionals are CMB random variables. We can use this fact to
* sample from CMM using a Gibbs sampler.
*/
// [[Rcpp::export("r_cmm_internal")]]
arma::mat r_cmm_internal_rcpp(unsigned int n, unsigned int m,
	const arma::vec& p, double nu, unsigned int burn, unsigned int thin,
	const arma::vec& x_init, unsigned int report);

//' @name cmm
//' @export
// [[Rcpp::export("d_cmm")]]
double d_cmm_rcpp(const arma::vec& x, const arma::vec& p, double nu,
	bool log = false, bool normalize = true);

//' @name cmm
//' @export
// [[Rcpp::export("normconst_cmm")]]
double normconst_cmm_rcpp(unsigned int m, const arma::vec& p, double nu,
	bool log = false);

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
// [[Rcpp::export("d_cmm_sample")]]
arma::vec d_cmm_sample_rcpp(const arma::mat& X, const arma::mat& P,
	const arma::vec& nu, bool log = false, bool normalize = true);

#endif
