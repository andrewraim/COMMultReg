#include "cmm.h"
#include "cmb.h"
#include "MultIterator.h"

/*
 * The CMM conditionals are CMB random variables. We can use this fact to
 * sample from CMM using a Gibbs sampler.
 */
arma::mat r_cmm_internal(unsigned int n, unsigned int m, const arma::vec& p,
	double nu, unsigned int burn, unsigned int thin,
	const arma::vec& x_init, unsigned int report_period)
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
		if ((r+1) % report_period == 0) {
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

double d_cmm(const arma::vec& x, const arma::vec& p, double nu, bool take_log,
	bool normalize)
{
	size_t k = x.n_elem;
	if (k != p.n_elem) {
		Rcpp::stop("Must have dim(x) == dim(p)");
	}

	size_t m = arma::sum(x);
	double ll = nu*lgamma(m+1) - nu*sum(lgamma(x+1)) + arma::dot(x, log(p));

	// Don't store the entire multinomial sample space; iterate through it.
	// Skip this if we don't need the probabilties to be normalized (e.g. if
	// it can be done more efficiently outside).
	if (normalize) {
		MultIterator itr(k, m);
		double sumprob = 0;
		for (; !itr.is_end(); itr.increment()) {
			// Note: tried this with ivec, but got different (wrong) answers
			const arma::vec& xx = itr.getCounts();
			sumprob += exp(nu*lgamma(m+1) - nu*sum(lgamma(xx+1)) + arma::dot(xx, log(p)));
		}
		ll -= log(sumprob);
	}

	if (take_log) { return ll; } else { return exp(ll); }
}

double normconst_cmm(unsigned int m, const arma::vec& p, double nu, bool take_log)
{
	size_t k = p.n_elem;
	arma::vec x = arma::zeros(k);
	x(0) = m;
	double ll = d_cmm(x, p, nu, true, false) - d_cmm(x, p, nu, true, true);
	if (take_log) { return ll; } else { return exp(ll); }
}

// A special vectorized version of the density that we use elsewhere in the package
arma::vec d_cmm_sample(const arma::mat& X, const arma::mat& P,
	const arma::vec& nu, bool take_log, bool normalize)
{
	unsigned int n = X.n_rows;
	unsigned int k = X.n_cols;
	if (n != P.n_rows || k != P.n_cols || n != nu.n_elem) {
		Rcpp::stop("Must have length(nu) == nrow(X) and dim(X) == dim(P)");
	}

	arma::vec out(n);
	for (size_t i = 0; i < n; i++) {
		out(i) = d_cmm(X.row(i).t(), P.row(i).t(), nu(i), take_log, normalize);
	}

	return out;
}
