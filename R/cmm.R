#' Conway Maxwell Multinomial distribution
#' 
#' Functions for the \eqn{\text{CMM}_k(m, \bm{p}, \nu)} distribution.
#' 
#' @param x A \eqn{k}-dimensional vector representing the outcome.
#' @param n Number of draws to produce.
#' @param m Number of trials in the CMB cluster.
#' @param p Probability parameter; a vector of positive numbers which sums to 1.
#' @param nu Dispersion parameter.
#' @param log `TRUE` or `FALSE`; if `TRUE`, return the value on the log-scale.
#' @param normalize `TRUE` or `FALSE`; if `FALSE`, do not
#' compute or apply the normalizing constant to each density value.
#' @param burn Number of initial draws to burn for Gibbs sampler.
#' @param thin Thinning interval for Gibbs sampler. A value of `s` means
#' that `s` iterations of the sampler will be carried out before saving each
#' of the `n` requested draws.
#' @param report How often to output progress for Gibbs sampler. A value of
#' `s` means that a progress message will be printed every `s` iterations
#' of the sampler.
#' @param x_init Initial value for Gibbs sampler. If `NULL`, it is set to
#' the extreme point \eqn{m \bm{e}_{\ell}} where \eqn{\ell = \text{argmax}_{j}\{p_j\}}. 
#' 
#' @return The values returned by each function are:
#' \itemize{
#' \item `d_cmm`: the CMM density \eqn{f(\bm{x} \mid m, \bm{p}, \nu)},
#' where \eqn{m} is assumed to be `sum(x)`.
#' \item `r_cmm`: an \eqn{n \times k} matrix of draws.
#' \item `normconst_cmm`: a number representing the normalizing constant \eqn{C(m, \bm{p}, \nu)}.
#' \item `e_cmm`: a \eqn{k}-dimensional vector representing \eqn{\text{E}(\bm{X})}.
#' \item `v_cmm`: a \eqn{k \times k} matrix representing \eqn{\text{Var}(\bm{X})}
#' }
#'
#' @details
#' Let \eqn{\Omega_{m,k}} denote the multinomial sample space based on \eqn{m} trials
#' and \eqn{k} categories. A random variable
#' \eqn{\bm{X} \sim \text{CMM}_k(m, \bm{p}, \nu)} has probability mass function
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
#' set.seed(1234)
#' m = 10
#' p = c(0.1,0.1,0.8)
#' nu = 0.8
#' 
#' x = r_cmm(100, m, p, nu, burn = 1000, thin = 10)
#' d_cmm(x[1,], p, nu, log = TRUE)
#' normconst_cmm(m, p, nu, log = TRUE)
#' e_cmm(m, p, nu)
#' v_cmm(m, p, nu)
#' 
#' @name cmm
NULL

#' @name cmm
#' @export
r_cmm = function(n, m, p, nu, burn = 0, thin = 1, x_init = NULL, report = NULL)
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
	r_cmm_internal(n, m, p, nu, burn, thin, x_init, report)
}

# Expected value of CMM, from the definition of E(X)
#' @name cmm
#' @export
e_cmm = function(m, p, nu)
{
	stopifnot(length(m) == 1)
	stopifnot(length(nu) == 1)
	stopifnot(all(p > 0))
	p = p / sum(p)

	k = length(p)
	xx = t(xsimplex(k, m))
	P_mat = matrix(p, nrow(xx), k, byrow = TRUE)
	nu_mat = rep(nu, nrow(xx))

	f_all_unnorm = d_cmm_sample(xx, P_mat, nu_mat, log = FALSE, normalize = FALSE)
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
	stopifnot(all(p > 0))
	p = p / sum(p)

	k = length(p)
	xx = t(xsimplex(k, m))
	P_mat = matrix(p, nrow(xx), k, byrow = TRUE)
	nu_mat = rep(nu, nrow(xx))

	f_all_unnorm = d_cmm_sample(xx, P_mat, nu_mat, log = FALSE, normalize = FALSE)
	f_all = as.numeric(normalize(f_all_unnorm))
	e = as.numeric(t(xx) %*% f_all)
	t(xx) %*% (f_all * xx) - e %*% t(e)
}
