#include <RcppArmadillo.h>
#include "COMMultReg.h"

arma::mat r_cmm_internal_rcpp(unsigned int n, unsigned int m,
	const arma::vec& p, double nu, unsigned int burn, unsigned int thin,
	const arma::vec& x_init, unsigned int report)
{
	return COMMultReg::r_cmm_internal(n, m, p, nu, burn, thin, x_init, report);
}

double d_cmm_rcpp(const arma::vec& x, const arma::vec& p, double nu, bool log,
	bool normalize)
{
	return COMMultReg::d_cmm(x, p, nu, log, normalize);
}

double normconst_cmm_rcpp(unsigned int m, const arma::vec& p, double nu, bool log)
{
	return COMMultReg::normconst_cmm(m, p, nu, log);
}

arma::vec d_cmm_sample_rcpp(const arma::mat& X, const arma::mat& P,
	const arma::vec& nu, bool log, bool normalize)
{
	return COMMultReg::d_cmm_sample(X,  P, nu, log, normalize);
}
