// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <map>
#include "MultIterator.h"

struct cmp_uvec {
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
typedef std::map<arma::uvec, unsigned int, cmp_uvec> freq_map_type;

//' @name cmm
//' @export
// [[Rcpp::export]]
Rcpp::List gunterize(const arma::umat& X, bool all = false)
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

