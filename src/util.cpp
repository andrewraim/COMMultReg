#include "util.h"

arma::mat pinv(const arma::mat& x) {
	return arma::pinv(x);
}

double logchoose(const arma::vec& x) {
	return lgamma(sum(x)+1) - sum(lgamma(x+1));
}

arma::vec normalize(const arma::vec& x, bool na_rm)
{
	size_t neg = 0;
	size_t pos = 0;
	size_t nan = 0;
	double ss = 0;
	for (size_t i = 0; i < x.n_elem; i++) {
		neg += (x(i) < 0);
		pos += (x(i) > 0);
		if (arma::is_finite(x(i))) {
			ss += x(i);
		} else {
			nan++;
		}
	}
	if (neg > 0) {
		throw std::range_error("x must have all non-negative components");
	}
	if (pos == 0) {
		throw std::range_error("x must have at least one positive component");
	}
	if (nan > 0 && !na_rm) {
		throw std::range_error("x has at least one NaN component which we cannot ignore");
	}

	return x / ss;
}

