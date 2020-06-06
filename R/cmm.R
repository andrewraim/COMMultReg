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
#' @param report_period 
#' @param x_init 
#' 
#' @return
#'
#' @examples
#' @name cmm
NULL

#' @name cmm
#' @export
r_cmm = function(n, m, p, nu, burn = 0, thin = 1, x_init = NULL, report_period = NULL)
{
	k = length(p)
	reps = burn + n*thin
	if (is.null(report_period)) {
		report_period = reps + 1
	}
	if (is.null(x_init)) {
		x_init = rep(0, k)
		idx_max = which.max(p)
		x_init[idx_max] = m
	}
	r_cmm_internal(n, m, p, nu, burn, thin, x_init, report_period)
}

# Expected value of CMM, from the definition of E(X)
#' @name cmm
#' @export
e_cmm = function(m, p, nu)
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
v_cmm = function(m, p, nu)
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
