#ifndef COMMULTREG_CMB_H
#define COMMULTREG_CMB_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace COMMultReg {

inline double d_cmb(unsigned int x, unsigned int m, double p,
	double nu, bool log = false, bool normalize = true)
{
	double logfx = nu*lgamma(m+1) - nu*lgamma(x+1) - nu*lgamma(m-x+1) +
		x*std::log(p) + (m-x)*std::log(1-p);

	if (normalize) {
		arma::vec z = arma::linspace<arma::vec>(0.0, m, m+1);
		arma::vec fz = exp(nu*lgamma(m+1) - nu*lgamma(z+1) - nu*lgamma(m-z+1) +
			z*std::log(p) + (m-z)*std::log(1-p));
		logfx -= std::log(sum(fz));
	}

	return log ? logfx : exp(logfx);	
}

inline arma::vec r_cmb(unsigned int n, unsigned int m, double p, double nu)
{
	arma::vec u = Rcpp::runif(n, 0, 1);
	arma::vec x(n);

	for (size_t i = 0; i < n; i++) {
		arma::vec z = arma::linspace<arma::vec>(0, m, m+1);
		arma::vec fz(m+1);
		for (unsigned int l = 0; l < m+1; l++) {
			fz(l) = d_cmb(z(l), m, p, nu, false, false);
		}
		fz = fz / arma::sum(fz);

		// This is a sneaky way of finding the max index such that:
		// u(i) > arma::cumsum(fz)
		// Note that it is possible for the RNG to draw exactly 1;
		// we have to handle this case specially.
		x(i) = (u(i) < 1)*sum(u(i) > arma::cumsum(fz)) + (u(i) >= 1)*m;
	}

	return x;	
}

inline double p_cmb(unsigned int x, unsigned int m, double p, double nu)
{
	double f_num = 0;
	double f_denom = 0;

	for (unsigned int i = 0; i < x+1; i++) {
		double fx = exp(nu*lgamma(m+1) - nu*lgamma(i+1) - nu*lgamma(m-i+1) +
			i*log(p) + (m-i)*log(1-p));
		f_num += fx;
		f_denom += fx;
	}

	for (unsigned int i = x+1; i < m+1; i++) {
		double fx = exp(nu*lgamma(m+1) - nu*lgamma(i+1) - nu*lgamma(m-i+1) +
			i*log(p) + (m-i)*log(1-p));
		f_denom += fx;
	}

	return f_num / f_denom;	
}

inline double q_cmb(unsigned int q, unsigned int m, double p, double nu)
{
	Rcpp::stop("This needs to be implemented");
	return -1;
}

inline double normconst_cmb(unsigned int m, double p, double nu,
	bool log = false)
{
	double log_f0 = d_cmb(0, m, p,nu, true, true);
	double out = -log_f0 + m*std::log(1-p);
	return log ? out : exp(out);
}

inline arma::vec d_cmb_sample(const arma::vec& x, const arma::vec& m,
	const arma::vec& p, const arma::vec& nu, bool log = false)
{
	unsigned int n = x.n_elem;
	arma::vec out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = d_cmb(x(i), m(i), p(i), nu(i), log);
	}

	return out;
}

}

#endif
