#include "cmb.h"

arma::vec r_cmb(unsigned int n, unsigned int m, double p, double nu)
{
	arma::vec u = Rcpp::runif(n, 0, 1);
	arma::vec x(n);
	
	for (size_t i = 0; i < n; i++) {
		arma::vec z = arma::linspace<arma::vec>(0, m, m+1);
		arma::vec fz(m+1);
		for (unsigned int l = 0; l < m+1; l++) {
			fz(l) = d_cmb(z(l), m, p, nu, false, false);
		}
		fz = fz / sum(fz);

		// This is a sneaky way of finding the max index such that:
		// u(i) > arma::cumsum(fz)
		// It is possible for the RNG to draw exactly 1. We have to handle this
		// case specially.
		x(i) = (u(i) < 1)*sum(u(i) > arma::cumsum(fz)) + (u(i) >= 1)*m;
	}

	return x;
}

double d_cmb(unsigned int x, unsigned int m, double p,
	double nu, bool take_log, bool normalize)
{
	double ll = nu*lgamma(m+1) - nu*lgamma(x+1) - nu*lgamma(m-x+1) +
		x*log(p) + (m-x)*log(1-p);

	if (normalize) {
		arma::vec z = arma::linspace<arma::vec>(0.0, m, m+1);
		arma::vec fz = exp(nu*lgamma(m+1) - nu*lgamma(z+1) - nu*lgamma(m-z+1) +
			z*log(p) + (m-z)*log(1-p));
		ll -= log(sum(fz));
	}

	if (take_log) { return ll; } else { return exp(ll); }
}

double d_cmb_normconst(unsigned int m, double p, double nu, bool take_log)
{
	double log_f0 = d_cmb(0, m, p,nu, true, true);
	double logC = -log_f0 + m*log(1-p);
	if (take_log) { return logC; } else { return exp(logC); }
}

double loglik_cmb(const arma::vec& y, const arma::vec& m,
	const arma::vec& p, const arma::vec& nu)
{
	unsigned int n = y.n_elem;
	double ll = 0;
	for (unsigned int i = 0; i < n; i++) {
		ll += d_cmb(y(i), m(i), p(i), nu(i), true);
	}

	return ll;
}

