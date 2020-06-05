#ifndef MULT_ITERATOR_H
#define MULT_ITERATOR_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

class MultIterator
{
public:
	MultIterator(unsigned int k, unsigned int m);
	virtual ~MultIterator();

	void increment();
	arma::uvec getState() const;
	arma::vec getCounts() const;
	bool is_end() const;

private:
	arma::uvec _x;
	arma::vec _xCounts;
	unsigned int _k;
	bool _end;
};

#endif

