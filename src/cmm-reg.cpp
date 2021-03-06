#include <map>
#include "cmm-reg.h"
#include "MultIterator.h"
#include "util.h"

// A comparator which specifies an ordering between two arma::uvec objects
struct comparator_uvec {
    bool operator()(const arma::uvec& a, const arma::uvec& b) const {
    	if (a.n_elem < b.n_elem) {
    		return true;
    	} else if (b.n_elem < a.n_elem) {
    		return false;
    	}

    	for (unsigned int i = 0; i < a.n_elem; i++) {
	    	if (a(i) < b(i)) {
	    		return true;
	    	} else if (b(i) < a(i)) {
	    		return false;
	    	}
    	}

        return false;
    }
};
typedef std::map<arma::uvec, unsigned int, comparator_uvec> freq_map_type;

Rcpp::List gunterize(const arma::umat& X, bool all)
{
	const arma::uvec& s = arma::vectorise(sum(X,1));
	unsigned int m = s(0);
	unsigned int k = X.n_cols;
	if (!arma::all(s == m)) {
		Rcpp::stop("All rows of X must have a common sum");
	}

	freq_map_type freq_map;
	if (all) {
		MultIterator itr(k, m);
		for (; !itr.is_end(); itr.increment()) {
			const arma::vec& v = itr.getCounts();
			const arma::uvec& u = arma::conv_to<arma::uvec>::from(v);
			freq_map[u] = 0;
		}
	}

	for (unsigned int i = 0; i < X.n_rows; i++) {
		freq_map[X.row(i).t()]++;
	}

	arma::umat Xnew(freq_map.size(), X.n_cols);
	arma::uvec freqs(freq_map.size());
	
	freq_map_type::const_iterator itr = freq_map.begin();
	for (unsigned int i = 0; itr != freq_map.end(); itr++, i++) {
		Xnew.row(i) = arma::trans(itr->first);
		freqs(i) = itr->second;
	}

	return Rcpp::List::create(
		Rcpp::Named("table") = Xnew,
		Rcpp::Named("freqs") = freqs
	);
}

Rcpp::List loglik_score_fim_cmm(const Rcpp::List& par,
	const Rcpp::List& dat_xform, unsigned int baseline)
{
	unsigned int L = dat_xform.length();
	const arma::vec& mu = par("mu");
	const arma::vec& gamma = par("gamma");
	const arma::vec& beta = par("beta");
	const arma::vec& vartheta = arma::join_cols(arma::join_cols(mu, gamma), beta);
	unsigned int qq = vartheta.n_elem;

	double ll = 0;
	arma::vec score = arma::zeros(qq, 1);
	arma::mat fim = arma::zeros(qq, qq);

	for (unsigned int l = 0; l < L; l++) {
		const Rcpp::List& dat = dat_xform(l);
		arma::rowvec ee = arma::zeros(1,L);
		ee(l) = 1;

		const arma::mat& zz = dat("z_obs");
		arma::mat zz_based = zz;
		zz_based.shed_col(baseline-1);
		unsigned int m = dat("m");
		const arma::vec& x = dat("x");
		const arma::vec& w = dat("w");
		const arma::vec& freqs = dat("freqs");
		unsigned int k = zz.n_cols;
		unsigned int n_l = sum(freqs);

		// This part of the calculation only requires the observed data
		for (unsigned int i = 0; i < zz.n_rows; i++) {
			arma::vec s = arma::join_rows(
				arma::join_rows(ee, logchoose(zz.row(i).t()) * w.t()),
				arma::kron(zz_based.row(i), x.t())
			).t();
			double svartheta = dot(s, vartheta) + log(n_l);

			ll += freqs(i) * svartheta - lgamma(freqs(i)+1);
			score = score + freqs(i) * s;
		}

		// This part of the calculation requires the whole mult sample space
		MultIterator itr(k, m);
		for (; !itr.is_end(); itr.increment()) {
			const arma::rowvec& z = itr.getCounts().t();
			arma::rowvec z_based = z;
			z_based.shed_col(baseline-1);

			arma::vec s = arma::join_rows(
				arma::join_rows(ee, logchoose(z.t()) * w.t()),
				arma::kron(z_based, x.t())
			).t();
			double svartheta = dot(s, vartheta) + log(n_l);
			double lambda = exp(svartheta);

			ll = ll - lambda;
			score = score - lambda * s;
			fim = fim + lambda * s * s.t();
		}

		// Check for user interrupt, since this function can take a long time.
		R_CheckUserInterrupt();
	}

	return Rcpp::List::create(
		Rcpp::Named("loglik") = ll,
		Rcpp::Named("score") = score,
		Rcpp::Named("fim") = fim);
}
