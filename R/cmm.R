#' Conway Maxwell Multinomial distribution
#' 
#' Functions for CMM distribution.
#' 
#' @param x 
#' @param m 
#' @param lambda 
#' @param nu 
#' @param log
#' @param take_log
#' @param normalize 
#' @param phi 
#' @param burn
#' @param R 
#' @param thin
#' @param report.period 
#' @param x.init 
#' 
#' @return
#'
#' @examples
#' @name cmm
NULL

# Expected value of CMM, from the definition of E(X)
#' @name cmm
#' @export
ecmm = function(m, p, nu)
{
	stopifnot(length(m) == 1)
	stopifnot(length(nu) == 1)
	stopifnot(is.integer(m))
	stopifnot(all(p > 0))
	p = p / sum(p)

	k = length(p)
	xx = t(xsimplex(k, m))
	p_mat = matrix(p, nrow(xx), k, byrow = TRUE)
	nu_mat = rep(nu, nrow(xx))
	m_mat = rep(m, nrow(xx))

	f_all_unnorm = d_cmm_vectorized(xx, m_mat, p_mat, nu_mat, take_log = FALSE, normalize = FALSE)
	f_all = normalize(f_all_unnorm)
	as.numeric(t(xx) %*% f_all)
}

# Variance of CMM, from the definition of Var(X)
#' @name cmm
#' @export
vcmm = function(m, p, nu)
{
	stopifnot(length(m) == 1)
	stopifnot(length(nu) == 1)
	stopifnot(is.integer(m))
	stopifnot(all(p > 0))
	p = p / sum(p)

	k = length(p)
	xx = t(xsimplex(k, m))
	p_mat = matrix(p, nrow(xx), k, byrow = TRUE)
	nu_mat = rep(nu, nrow(xx))
	m_mat = rep(m, nrow(xx))

	f_all_unnorm = d_cmm_vectorized(xx, m_mat, p_mat, nu_mat, take_log = FALSE, normalize = FALSE)
	f_all = as.numeric(normalize(f_all_unnorm))
	e = as.numeric(t(xx) %*% f_all)
	t(xx) %*% (f_all * xx) - e %*% t(e)
}
