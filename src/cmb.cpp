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
		// Note that it is possible for the RNG to draw exactly 1;
		// we have to handle this case specially.
		x(i) = (u(i) < 1)*sum(u(i) > arma::cumsum(fz)) + (u(i) >= 1)*m;
	}

	return x;
}

double d_cmb(unsigned int x, unsigned int m, double p,
	double nu, bool take_log, bool normalize)
{
	double logfx = nu*lgamma(m+1) - nu*lgamma(x+1) - nu*lgamma(m-x+1) +
		x*log(p) + (m-x)*log(1-p);

	if (normalize) {
		arma::vec z = arma::linspace<arma::vec>(0.0, m, m+1);
		arma::vec fz = exp(nu*lgamma(m+1) - nu*lgamma(z+1) - nu*lgamma(m-z+1) +
			z*log(p) + (m-z)*log(1-p));
		logfx -= log(sum(fz));
	}

	if (take_log) { return logfx; } else { return exp(logfx); }
}

double normconst_cmb(unsigned int m, double p, double nu, bool take_log)
{
	double log_f0 = d_cmb(0, m, p,nu, true, true);
	double logC = -log_f0 + m*log(1-p);
	if (take_log) { return logC; } else { return exp(logC); }
}

arma::vec d_cmb_sample(const arma::vec& x, const arma::vec& m,
	const arma::vec& p, const arma::vec& nu, bool take_log)
{
	unsigned int n = x.n_elem;
	arma::vec out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = d_cmb(x(i), m(i), p(i), nu(i), take_log);
	}

	return out;
}

