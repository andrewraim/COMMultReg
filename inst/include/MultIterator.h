#ifndef MULT_ITERATOR_H
#define MULT_ITERATOR_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace COMMultReg {

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

inline MultIterator::MultIterator(unsigned int k, unsigned int m)
	: _x(m), _xCounts(k), _k(k), _end(false)
{
	if (m == 0) {
		Rcpp::stop("k must be at least 1");
	}

	// Initialize all trials to be in the 0 position.
	// Accordingly, the counts should be (m, 0, ..., 0),
	_x.fill(0);
	_xCounts.fill(0);
	_xCounts(0) = m;
}

inline MultIterator::~MultIterator()
{
}

/*
* The logic here is from the function gsl_multiset_next in the
* GSL package; see the file multiset.c in its source code.
* 
* Be careful - this can get confusing! A multiset describes the
* positions of the trials. A multinomial outcome counts the number
* of trials in each position. 
* 
* We try to update the multinomial outcome at the same time as the
* counts here, even though we could convert the multiset to a
* multinomial outcome later. We also save to an arma::vec rather
* than an arma::vec. Both of these are done specifically to speed up
* our application.
*/
inline void MultIterator::increment()
{
	unsigned int m = _x.size();
	size_t i = m - 1;

	while (i > 0 && _x(i) == _k-1) {
		--i;
	}

	if (i == 0 && _x(0) == _k-1) {
		_end = true;
	} else {
		_xCounts(_x(i))--;
		_x(i)++;
		_xCounts(_x(i))++;
		while (i < m - 1) {
			_xCounts(_x(i+1))--;
			_x(i+1) = _x(i);
			_xCounts(_x(i+1))++;
			++i;
		}
	}
}

inline arma::uvec MultIterator::getState() const
{
	return _x;
}

inline arma::vec MultIterator::getCounts() const
{
	return _xCounts;
}

inline bool MultIterator::is_end() const
{
	return _end;
}

}

#endif
