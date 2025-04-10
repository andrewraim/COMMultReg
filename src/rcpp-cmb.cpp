#include "COMMultReg.h"

double d_cmb_rcpp(unsigned int x, unsigned int m, double p, double nu,
	bool log, bool normalize)
{
	return COMMultReg::d_cmb(x, m, p, nu, log, normalize);
}

arma::vec r_cmb_rcpp(unsigned int n, unsigned int m, double p, double nu)
{
	return COMMultReg::r_cmb(n, m, p, nu);
}

double p_cmb_rcpp(unsigned int x, unsigned int m, double p, double nu)
{
	return COMMultReg::p_cmb(x, m, p, nu);
}

double q_cmb_rcpp(unsigned int q, unsigned int m, double p, double nu)
{
	return COMMultReg::q_cmb(q, m, p, nu);
}

double normconst_cmb_rcpp(unsigned int m, double p, double nu, bool log)
{
	return COMMultReg::normconst_cmb(m, p, nu, log);
}

arma::vec d_cmb_sample_rcpp(const arma::vec& x, const arma::vec& m,
	const arma::vec& p, const arma::vec& nu, bool log)
{
	return COMMultReg::d_cmb_sample(x, m, p, nu, log);
}
