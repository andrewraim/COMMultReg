#' Conway Maxwell Binomial distribution
#' 
#' Functions for CMB distribution.
#' 
#' @param x 
#' @param m 
#' @param p 
#' @param nu 
#' @param take_log
#' @param normalize 
#' @param phi 
#' 
#' @return
#'
#' @examples
#' @name cmb
NULL

#' @name cmb
#' @export
dcmb.old <- function(x, m, p, nu, take_log = FALSE, normalize = TRUE)
{
	stopifnot(all(is.integer(x)))
	stopifnot(all(is.integer(m)))
	stopifnot(all(0 < p & p < 1))

	n <- length(x)
	if (length(m) == 1) { m <- rep(m, n) }
	if (length(p) == 1) { p <- rep(p, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }
	fx <- dcmb_cpp(as.numeric(x), as.numeric(m), p, nu, take_log, normalize)
	as.numeric(fx)
}

#' @name cmb
#' @export
rcmb.old <- function(n, m, p, nu)
{
	stopifnot(is.integer(m))
	stopifnot(all(is.integer(m)))
	stopifnot(all(0 < p & p < 1))

	if (length(m) == 1) { m <- rep(m, n) }
	if (length(p) == 1) { p <- rep(p, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }
	x <- rcmb_cpp(n, as.numeric(m), p, nu)
	as.integer(x)
}

# A slow / pure-R version of CMB density, purely for comparison
#' @name cmb
#' @export
dcmb.slow <- function(x, m, p, nu, take_log = FALSE, normalize = TRUE)
{
	n <- length(x)
	if (length(m) == 1) { m <- rep(m, n) }
	if (length(p) == 1) { p <- rep(p, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }
	exp(nu*lgamma(m+1) - nu*lgamma(x+1) - nu*lgamma(m-x+1) +
		x*log(p) + (m-x)*log(1-p))
}

# Expected value of CMB, from the definition of E(X)
#' @name cmb
#' @export
ecmb <- function(m, p, nu)
{
	stopifnot(length(m) == 1)
	stopifnot(length(p) == 1)
	stopifnot(length(nu) == 1)
	stopifnot(is.integer(m))

	xx <- 0:m
	f.all.unnorm <- numeric(m+1)
	for (i in seq_along(xx)) {
		f.all.unnorm[i] <- d_cmb(xx[i], m, p, nu, take_log = FALSE, normalize = FALSE)
	}

	f.all <- normalize(f.all.unnorm)
	sum(xx * f.all)
}

# Expected 2nd factorial moment of CMB, from the definition of E[X(X-1)]
#' @name cmb
#' @export
ecmb_factorial2 <- function(m, p, nu)
{
	stopifnot(length(m) == 1)
	stopifnot(length(p) == 1)
	stopifnot(length(nu) == 1)
	stopifnot(is.integer(m))

	xx <- 0:m
	f.all.unnorm <- numeric(m+1)
	for (i in seq_along(xx)) {
		f.all.unnorm[i] <- d_cmb(xx[i], m, p, nu, take_log = FALSE, normalize = FALSE)
	}

	f.all <- normalize(f.all.unnorm)
	sum(xx * (xx-1) * f.all)
}

