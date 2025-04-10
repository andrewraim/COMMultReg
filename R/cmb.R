#' Conway Maxwell Binomial distribution
#' 
#' Functions for \eqn{\text{CMB}(m, p, \nu)} distribution.
#' 
#' @param x A scalar representing the outcome.
#' @param n Number of draws to produce.
#' @param m Number of trials in the CMB cluster.
#' @param p Probability parameter
#' @param nu Dispersion parameter
#' @param log `TRUE` or `FALSE`; if `TRUE`, return
#' the value on the log-scale.
#' @param normalize `TRUE` or `FALSE`; if `FALSE`, do not
#' compute or apply the normalizing constant to each density value.
#' 
#' @return The values returned by each function are:
#' \itemize{
#' \item `d_cmb`: a number representing the CMB density \eqn{f(x \mid m, p, \nu)}.
#' \item `r_cmb`: an \eqn{n}-dimensional vector of draws.
#' \item `normconst_cmb`: a number representing the normalizing constant \eqn{C(m, p, \nu)}.
#' \item `e_cmb`: a number representing \eqn{\text{E}(X)}.
#' \item `v_cmb`: a number representing \eqn{\text{Var}(X)}
#' }
#' 
#' @details
#' A random variable
#' \eqn{X \sim \text{CMB}(m, p, \nu)} has probability mass function
#' \deqn{
#' f(x \mid m, p, \nu) = C(m, p, \nu)^{-1} {m \choose x}^\nu p^{x} (1-p)^{m-x},
#' \quad x \in \{0, 1, \ldots, m\}
#' }
#' with
#' \deqn{
#' C(m, p, \nu) = \sum_{x = 0}^m
#' {m \choose x}^\nu p^{x} (1-p)^{m-x}
#' }
#' as the normalizing constant.
#' 
#' CMB can be considered a two-dimensional case of CMM. Furthermore, CMB is
#' useful in drawing from the general CMM distribution; see
#' Morris, Raim, and Sellers (2020+).
#' 
#' @examples
#' set.seed(1234)
#' m = 10
#' p = 0.7
#' nu = 0.8
#' 
#' x = r_cmb(100, m, p, nu)
#' d_cmb(x[1], m, p, nu, log = TRUE)
#' normconst_cmb(m, p, nu, log = TRUE)
#' e_cmb(m, p, nu)
#' v_cmb(m, p, nu)
#' 
#' @name cmb
NULL

# Expected value of CMB, computed from the definition of E(X)
#' @name cmb
#' @export
e_cmb = function(m, p, nu)
{
	stopifnot(length(m) == 1)
	stopifnot(length(p) == 1)
	stopifnot(length(nu) == 1)

	xx = 0:m
	f_all_unnorm = numeric(m+1)
	for (i in seq_along(xx)) {
		f_all_unnorm[i] = d_cmb(xx[i], m, p, nu, log = FALSE, normalize = FALSE)
	}

	f_all = normalize(f_all_unnorm)
	sum(xx * f_all)
}

# Variance of CMB, computed from the definition of Var(X)
#' @name cmb
#' @export
v_cmb = function(m, p, nu)
{
	stopifnot(length(m) == 1)
	stopifnot(length(p) == 1)
	stopifnot(length(nu) == 1)

	xx = 0:m
	f_all_unnorm = numeric(m+1)
	for (i in seq_along(xx)) {
		f_all_unnorm[i] = d_cmb(xx[i], m, p, nu, log = FALSE, normalize = FALSE)
	}

	f_all = normalize(f_all_unnorm)
	sum(xx^2 * f_all) - sum(xx * f_all)^2
}
