#ifndef CMB_H
#define CMB_H

//' @name cmb
//' @export
// [[Rcpp::export("d_cmb")]]
double d_cmb_rcpp(unsigned int x, unsigned int m, double p,
	double nu, bool log = false, bool normalize = true);

//' @name cmb
//' @export
// [[Rcpp::export("r_cmb")]]
arma::vec r_cmb_rcpp(unsigned int n, unsigned int m, double p, double nu);

//' @name cmb
//' @export
// [[Rcpp::export("p_cmb")]]
double p_cmb_rcpp(unsigned int x, unsigned int m, double p, double nu);

//' @name cmb
//' @export
// [[Rcpp::export("q_cmb")]]
double q_cmb_rcpp(unsigned int q, unsigned int m, double p, double nu);

//' @name cmb
//' @export
// [[Rcpp::export("normconst_cmb")]]
double normconst_cmb_rcpp(unsigned int m, double p, double nu, bool log = false);

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
//' @param log \code{TRUE} or \code{FALSE}; if \code{TRUE}, return the
//' value on the log-scale.
//' 
//' @return
//' A vector of density values
//' \eqn{f(x_1 \mid m_1, p_1, \nu_1), \ldots, f(x_n \mid m_n, p_n, \nu_n),}
//' which may be on the log-scale and/or unnormalized
//' according to input arguments. See \link{cmb}.
//'
//' @examples
//' set.seed(1234)
//' 
//' n = 20
//' m = rep(10, n)
//' 
//' x = rnorm(n)
//' X = model.matrix(~ x)
//' beta = c(-1, 1)
//' p = plogis(X %*% beta)
//' 
//' w = rnorm(n)
//' W = model.matrix(~ w)
//' gamma = c(0.1, -0.1)
//' nu = X %*% gamma
//' 
//' y = numeric(n)
//' for (i in 1:n) {
//'     y[i] = r_cmb(1, m[i], p[i], nu[i])
//' }
//' 
//' d_cmb_sample(y, m, p, nu, log = TRUE)
//' 
//' @export
// [[Rcpp::export("d_cmb_sample")]]
arma::vec d_cmb_sample_rcpp(const arma::vec& x, const arma::vec& m,
	const arma::vec& p, const arma::vec& nu, bool log = false);

#endif
