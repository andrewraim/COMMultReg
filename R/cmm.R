#' Conway Maxwell Multinomial distribution
#' 
#' Functions for the \eqn{\textrm{CMM}_k(m, \bm{p}, \nu)} distribution.
#' 
#' @param x A eqn{k}-dimensional vector representing the outcome.
#' @param n Number of draws to produce.
#' @param m Number of trials in the CMB cluster.
#' @param p Probability parameter; a vector of positive numbers which sums to 1.
#' @param nu Dispersion parameter.
#' @param take_log \code{TRUE} or \code{FALSE}; if \code{TRUE}, return
#' the value on the log-scale.
#' @param normalize \code{TRUE} or \code{FALSE}; if \code{FALSE}, do not
#' compute or apply the normalizing constant to each density value.
#' @param burn Number of initial draws to burn for Gibbs sampler.
#' @param thin Thinning interval for Gibbs sampler. A value of \code{s} means
#' that \code{s} iterations of the sampler will be carried out before saving each
#' of the \code{n} requested draws.
#' @param report_period How often to output progress for Gibbs sampler. A value of
#' \code{s} means that a progress message will be printed every \code{s} iterations
#' of the sampler.
#' @param x_init TBD
#' 
#' @return The values returned by each function are:
#' \itemize{
#' \item \code{d_cmm}: a number representing the CMM density \eqn{f(\bm{x} \mid m, \bm{p}, \nu)}.
#' \item \code{r_cmm}: an \eqn{n \times k} matrix of draws.
#' \item \code{normconst_cmm}: a number representing the normalizing constant \eqn{C(m, \bm{p}, \nu)}.
#' \item \code{e_cmm}: a \eqn{k}-dimensional vector representing \eqn{\textrm{E}(\bm{X})}.
#' \item \code{v_cmm}: a \eqn{k \times k} matrix representing \eqn{\textrm{Var}(\bm{X})}
#' }
#'
#' @details
#' Let \eqn{\Omega_{m,k}} denote the multinomial sample space. A random variable
#' \eqn{\bm{X} \sim \textrm{CMM}_k(m, \bm{p}, \nu)} has probability mass function
#' \deqn{
#' f(\bm{x} \mid m, \bm{p}, \nu) = C(m, \bm{p}, \nu)^{-1} {m \choose x_1 \cdots x_k}^\nu
#' p_1^{x_1} \cdots p_k^{x_k}, \quad \bm{x} \in \Omega_{m,k}
#' }
#' with
#' \deqn{
#' C(m, \bm{p}, \nu) = \sum_{\bm{x} \in \Omega_{m,k}}
#' {m \choose x_1 \cdots x_k}^\nu
#' p_1^{x_1} \cdots p_k^{x_k}
#' }
#' as the normalizing constant.
#' 
#' This package computes CMM density values in C++ to improve the performance of iterating
#' over the sample space. A Gibbs sampler is used to draw samples from CMM. Further details
#' can be found in Morris, Raim, and Sellers (2020+).
#'
#' @examples
#' stop("TBD")
#' 
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

	f_all_unnorm = d_cmm_sample(xx, m_mat, p_mat, nu_mat, take_log = FALSE, normalize = FALSE)
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

	f_all_unnorm = d_cmm_sample(xx, m_mat, p_mat, nu_mat, take_log = FALSE, normalize = FALSE)
	f_all = as.numeric(normalize(f_all_unnorm))
	e = as.numeric(t(xx) %*% f_all)
	t(xx) %*% (f_all * xx) - e %*% t(e)
}
