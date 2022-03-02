#' Conway Maxwell Multinomial Regression
#' 
#' Functions for CMM Regression.
#' 
#' @param formula_x Regression formula to determine \eqn{\bm{X}} design matrix.
#' See details. The outcome must be specified here.
#' @param formula_w Regression formula to determine \eqn{\bm{W}} design matrix.
#' @param data An optional \code{data.frame} with variables to be used with
#' regression formulas. Variables not found here are read from the environment.
#' @param object An object output from \code{cmm_reg} or \code{cmm_reg_raw}.
#' @param x Used with the \code{print} function; same meaning as \code{object}
#' argument.
#' @param y An \eqn{n \times k} matrix of outcomes, where the \eqn{i}th row
#' \eqn{\bm{x}_i^\top} represents the \eqn{i}th observation.
#' @param X The \eqn{\bm{X}} design matrix; see details.
#' @param W The \eqn{\bm{W}} design matrix; see details.
#' @param extended boolean; if \code{FALSE}, drop terms associated with the
#' extended parameters \eqn{\mu}. See details.
#' @param beta_init A \eqn{p_x \times (k-1)} matrix whose columns correspond to
#' \eqn{\bm{\beta_1}, \ldots, \bm{\beta_{k-1}}}. See details. If provided,
#' will be used for the starting value in Newton-Raphson. If \code{NULL},
#' a default value will be used.
#' @param gamma_init A vector which serves the starting value for
#' \eqn{\gamma}, if provided. Otherwise, if \code{NULL}, a default value
#' will be used.
#' @param control a list of control parameters; see \link{control}.
#' @param k the penalty per parameter to use in AIC; default is 2.
#' @param ... Additional optional arguments
#' 
#' @return An object of class \code{cmm_reg} containing the result.
#' 
#' @details
#' Fit the model
#' \deqn{
#' \bm{Y}_i \sim \textrm{CMM}_k(m_i, \bm{p}_i, \nu_i),
#' \quad i = 1, \ldots, n
#' }
#' assuming multinomial logit link
#' \deqn{
#' \log\left(\frac{p_{i2}}{p_{i1}} \right) = \bm{x}_i^\top \bm{\beta}_1,
#' \; \ldots, \;
#' \log\left(\frac{p_{i,k} }{ p_{i1}} \right) = \bm{x}_i^\top \bm{\beta}_{k-1},
#' }
#' and identity link
#' \deqn{
#' \nu_i = \bm{w}_i^\top \bm{\gamma}.
#' }
#' The first category is assumed to be the baseline by default, but this can be
#' changed to category \code{b} by specifying the \code{base = b} argument.
#' 
#' The function \code{cmm_reg} provides a formula interface, while
#' \code{cmm_reg_raw} provides a "raw" interface where variables \code{y},
#' \code{X}, and \code{W} can be provided directly.
#' 
#' Fitting is carried out with a Newton-Raphson algorithm on the extended parameter
#' \eqn{\bm{\vartheta} = (\bm{\mu}, \bm{\psi})}, where
#' \eqn{\bm{\psi} = (\bm{\gamma}, \bm{\beta}_1, \ldots, \bm{\beta}_{k-1})} are
#' the regression coefficients of interest and \eqn{\bm{\mu}}
#' contains elements of the form \eqn{-\log C(\bm{p}, \nu; m) + m \log p_b} which are
#' typically not of direct interest to the analyst. See Morris, Raim, and
#' Sellers (2020+) for further details.
#' 
#' Let \eqn{\bm{\vartheta}^{(g)}} denote the estimate from the \eqn{g}th iteration.
#' The algorithm is considered to have converged when 
#' \eqn{\sum_{j} |\vartheta_j^{(g)} - \vartheta_j^{(g-1)}| < \epsilon}
#' is sufficiently small, or failed to converge when a maximum number of iterations
#' has been reached. These values can be specified via the \code{control} argument.
#'
#' Several auxiliary functions are provided for convenience:
#' \itemize{
#' \item \code{cmm_reg_control} provides a convenient way to construct the
#' \code{control} argument.
#' \item \code{logLik} returns the log-likelihood at the solution
#' \eqn{\hat{\bm{\psi}}}.
#' \item \code{AIC} returns the Akaike information criterion at the solution.
#' \item \code{coef} returns a list with elements \code{beta} and \code{gamma}
#' at the solution \eqn{\hat{\bm{\psi}}}. If \code{extended = TRUE}, an element
#' with \eqn{\mu} is also returned.
#' \item \code{vcov} returns an estimate of the covariance matrix of
#' \eqn{\hat{\bm{\psi}}} based on the information matrix. If
#' \code{extended = TRUE}, elements corresponding to \eqn{\mu} are also
#' included.
#' \item \code{print} displays estimates and standard errors.
#' }
#'
#' @examples
#' # Generate data from CMM regression model
#' set.seed(1234)
#' 
#' n = 200
#' m = rep(10, n)
#' k = 3
#' 
#' x = rnorm(n)
#' X = model.matrix(~ x)
#' beta = matrix(NA, 2, k-1)
#' beta[1,] = -1
#' beta[2,] = 1
#' P = t(apply(X %*% beta, 1, inv_mlogit))
#' 
#' w = rnorm(n)
#' W = model.matrix(~ w)
#' gamma = c(1, -0.1)
#' nu = X %*% gamma
#' 
#' y = matrix(NA, n, k)
#' for (i in 1:n) {
#'     y[i,] = r_cmm(1, m[i], P[i,], nu[i], burn = 200)
#' }
#' 
#' dat = data.frame(y1 = y[,1], y2 = y[,2], y3 = y[,3], x = x, w = w)
#' 
#' # Fit CMM regression with formula interface
#' cmm_out = cmm_reg(formula_x = y ~ x, formula_w = ~w, data = dat)
#' print(cmm_out)
#' logLik(cmm_out)
#' AIC(cmm_out)
#' coef(cmm_out)
#' vcov(cmm_out)
#' 
#' # Alteratively, use the "raw" interface
#' cmm_raw_out = cmm_reg_raw(y, X, W)
#' print(cmm_raw_out)
#' 
#' @name cmm_reg
NULL

#' @name cmm_reg
#' @export
cmm_reg = function(formula_x, formula_w = ~ 1, data = NULL,
	beta_init = NULL, gamma_init = NULL, control = NULL, ...)
{
	# Parse formula_x. The response should be specified here.
	mf = model.frame(formula_x, data, ...)
	y = model.response(mf)
	X = model.matrix(formula_x, mf)
	p_x = ncol(X)

	off_x = model.offset(mf)
	if(!is.null(off_x)) {
		stop("offset in formula is currently not supported")
	}
	weights = model.weights(mf)
	if(!is.null(weights)) {
		stop("weights argument is currently not supported")
	}

	# Parse formula_w
	mf = model.frame(formula_w, data, ...)
	W = model.matrix(formula_w, mf)
	p_w = ncol(W)

	off_w = model.offset(mf)
	if(!is.null(off_w)) {
		stop("offset in formula is currently not supported")
	}

	n = nrow(y)
	k = ncol(k)

	cmm_reg_raw(y = y, X = X, W = W, beta_init = beta_init,
		gamma_init = gamma_init, control = control)
}

#' @name cmm_reg
#' @export
cmm_reg_raw = function(y, X, W, beta_init = NULL, gamma_init = NULL, control = NULL)
{
	k = ncol(y)
	p_x = ncol(X)
	p_w = ncol(W)

	dat_xform = transform_data(y, X, W)
	L = length(dat_xform)

	if (is.null(control)) {
		control = cmm_reg_control()
	}
	if (is.null(beta_init)) {
		beta_init = matrix(0, p_x, k-1)
	}
	if (is.null(gamma_init)) {
		gamma_init = numeric(p_w)
	}

	# Set mu_init to something consistent with other initial values
	mu_init = extended_intercepts(dat_xform, control$base, beta_init, gamma_init)

	par_init = list(mu = mu_init, gamma = gamma_init, beta = beta_init)
	out = newton_raphson(par_init = par_init, dat_xform = dat_xform,
		max_iter = control$max_iter, verbose = control$verbose,
		base = control$base, tol = control$tol)
	
	# Make some adjustments to the Newton Raphson output and return it
	class(out) = "cmm_reg"
	out$loglik_poisson = out$loglik
	out$base = control$base
	out$loglik = loglik_cmm_xform(out$par, dat_xform, base = control$base)
	out$y = y
	out$X = X
	out$W = W

	# Make sure labels for results are set properly
	out$x_names = colnames(X)
	colnames(out$par$beta) = sprintf("cat%d", (1:k)[-control$base])
	if (is.null(out$x_names)) {
		rownames(out$par$beta) = sprintf("X%d", 1:p_x)
	} else if (length(out$x_names) != p_x) {
		warning("length of x_names does not match ncol(X)")
	} else {
		rownames(out$par$beta) = out$x_names
	}

	out$w_names = colnames(W)
	names(out$par$gamma) = out$w_names

	return(out)
}

#' Conway Maxwell Multinomial Regression Control
#' 
#' An object with control arguments for CMM regression.
#' 
#' @param base in an integer which specifies which category is considered
#' the baseline. The default is 1.
#' @param tol specifies the convergence tolerance \eqn{\epsilon}. The
#' default is \code{1e-8}. 
#' @param verbose is a boolean; if \code{TRUE}, print informative
#' messages during fitting. Default is \code{FALSE}.
#' @param max_iter specifies the maximum number of Newton-Raphson
#' iterations. Default is \code{200}.
#'
#' @export
cmm_reg_control = function(base = 1, tol = 1e-8, verbose = FALSE, max_iter = 200)
{
	list(base = base, tol = tol, verbose = verbose, max_iter = max_iter)
}

par2vec = function(par)
{
	as.numeric(c(par$mu, par$gamma, as.numeric(par$beta)))
}

vec2par = function(coefs, L, p_x, p_w, k)
{
	list(
		mu = coefs[1:L],
		gamma = coefs[L + 1:p_w],
		beta = matrix(coefs[L + p_w + seq_len(p_x*(k-1))], p_x, k-1)
	)
}

transform_data = function(y, X, W, tol = 1e-8)
{
	n = nrow(y)
	k = ncol(y)
	m = rowSums(y)
	p_x = ncol(X)
	p_w = ncol(W)
	if (length(m) == 1) { m = rep(m,n) }
	stopifnot(n == length(m))
	stopifnot(n == nrow(X))
	stopifnot(n == nrow(W))

	# There's probably a more efficient way to do this, but we shouldn't need to
	# repeat it very often.
	vals = cbind(m, X, W)
	vals_uniq = unique(vals)
	L = nrow(vals_uniq)
	dat_grp = list()
	for (l in 1:L) {
		diff = matrix(vals_uniq[l,], nrow = nrow(vals), ncol = ncol(vals), byrow = TRUE) - vals
		idx_match = which(sqrt(rowSums(diff^2)) < tol)
		idx1 = idx_match[1]

		dat_grp[[l]] = list()
		dat_grp[[l]]$idx = idx_match
		dat_grp[[l]]$m = m[idx1]
		dat_grp[[l]]$x = t(X[idx1,])
		dat_grp[[l]]$w = t(W[idx1,])

		if (length(idx_match) > 1) {
			y_match = y[idx_match,]
		} else {
			y_match = t(y[idx_match,])
		}

		# Partition the observations into their values and frequencies
		out = gunterize(y_match, all = FALSE)
		dat_grp[[l]]$z_obs = out$table
		dat_grp[[l]]$freqs = as.numeric(out$freqs)
	}

	return(dat_grp)
}

# Compute the extended intercepts ("mu") for a given beta and gamma.
# These depend on the data, which we take in transformed format.
extended_intercepts = function(dat_xform, base, beta, gamma)
{
	L = length(dat_xform)
	mu = numeric(L)

	for (l in 1:L) {
		dat = dat_xform[[l]]
		m = dat$m
		p = inv_mlogit(dat$x %*% beta)
		nu = dat$w %*% gamma
		mu[l] = -normconst_cmm(m, p, nu, take_log = TRUE) + m*log(p[base])
	}

	return(mu)
}

loglik_cmm_xform = function(par, dat_xform, base = 1)
{
	L = length(dat_xform)
	coefs = par2vec(par)
	qq = length(coefs)
	ll = 0

	# Compute the loglik from the Poisson regression form.
	# This should be much faster than computing the log-likelihood from scratch,
	# since the normalizing constants are in the parameters.
	for (l in 1:L) {
		dat = dat_xform[[l]]
		ee = numeric(L)
		ee[l] = 1
		n_l = sum(dat_xform[[l]]$freqs)

		for (i in 1:nrow(dat$z_obs)) {
			s = c(ee,
				logchoose(dat$z_obs[i,]) * t(dat$w),
				t(dat$z_obs[i,-base]) %x% t(dat$x)
			)
			scoefs = sum(s * coefs)
			ll = ll + dat$freqs[i] * scoefs
		}
	}

	return(ll)
}

newton_raphson = function(par_init, dat_xform, base = 1, tol = 1e-6,
	max_iter = 100, verbose = FALSE)
{
	L = length(dat_xform)
	p_x = nrow(par_init$beta)
	p_w = length(par_init$gamma)
	k = ncol(par_init$beta) + 1

	delta = Inf
	coefs = par2vec(par_init)
	iter = 0
	st = Sys.time()

	if (verbose) {
		out = loglik_score_fim_cmm(par_init, dat_xform, baseline = base)
		logger("iter 0: Poisson loglik: %g\n", out$loglik)
		print(coefs[-(1:L)])
	}

	while (delta > tol && iter < max_iter) {
		iter = iter + 1
		coefs_old = coefs

		par = vec2par(coefs, L, p_x, p_w, k)
		out = loglik_score_fim_cmm(par, dat_xform, baseline = base)
		move = tryCatch({
			solve(out$fim, out$score)
		}, error = function(e) {
			warning("FIM was close to singular, trying pseudo-inverse instead...")
			pinv(out$fim) %*% out$score
		})
		coefs = coefs + move
		delta = sum(abs(coefs - coefs_old))

		if (verbose) {
			logger("iter %d: Poisson loglik: %g delta: %g\n", iter, out$loglik, delta)
			print(coefs[-(1:L)])
		}
	}

	elapsed_sec = as.numeric(Sys.time() - st, units = "secs")

	ret = list(par = par, loglik = out$loglik, score = out$score,
		fim = out$fim, tol = delta, iter = iter, L = L, elapsed_sec = elapsed_sec)
	return(ret)
}

#' @name cmm_reg
#' @export
logLik.cmm_reg = function(object, ...)
{
	object$loglik
}

#' @name cmm_reg
#' @export
AIC.cmm_reg = function(object, ..., k = 2)
{
	p = length(object$par$gamma) + length(object$par$beta)
	-2*logLik(object) + k*p
}

#' @name cmm_reg
#' @export
coef.cmm_reg = function(object, extended = FALSE, ...)
{
	if (extended) {
		ret = object$par
	} else {
		ret = list(gamma = object$par$gamma, beta = object$par$beta)
	}
	return(ret)
}

#' @name cmm_reg
#' @export
vcov.cmm_reg = function(object, extended = FALSE, ...)
{
	coefs = par2vec(object$par)	
	idx = setdiff(1:length(coefs), 1:object$L)
	V = tryCatch({
		solve(object$fim)
	}, error = function(e) {
		warning("FIM was close to singular, trying pseudo-inverse instead...")
		pinv(object$fim)
	})

	par_names = get_par_names(object, extended = TRUE)
	rownames(V) = par_names
	colnames(V) = par_names

	if (extended) {
		return(V)	
	} else {	
		return(V[idx,idx])	
	}

	return(V)
}

#' @export
fitted.cmm_reg = function(object, newX, newW, ...)
{
	par = object$par
	k = ncol(par$beta) + 1
	p_x = nrow(par$beta)
	p_w = nrow(par$gamma)
	base = object$base

	stopifnot(ncol(newX) == p_x)
	stopifnot(ncol(newW) == p_w)

	dat = data.frame(
		t(apply(newX %*% par$beta, 1, inv_mlogit, base = base)),
		newW %*% par$gamma
	)
	colnames(dat) = c(sprintf("p%d", 1:k), "nu")
	return(dat)
}

get_par_names = function(object, extended = FALSE, ...)
{
	par = coef(object, extended = TRUE)
	L = length(par$mu)
	k = ncol(par$beta) + 1
	p_x = nrow(par$beta)
	p_w = length(par$gamma)

	x_names = object$x_names
	cats = setdiff(1:k, object$base)
	if (is.null(x_names)) {
		label_var = rep(1:p_x, times = k-1)
		label_dim = rep(cats, each = p_x)
		snames = sprintf("%d:beta%d", label_dim, label_var)
	} else {
		stopifnot(length(x_names) == nrow(par$beta))
		label_var = rep(x_names, times = k-1)
		label_dim = rep(cats, each = p_x)
		snames = sprintf("%d:%s", label_dim, x_names)
	}

	w_names = object$w_names
	if (is.null(w_names)) {
		w_names = sprintf("gamma%d", 1:length(par$gamma))
	} else {
		stopifnot(length(w_names) == length(par$gamma))
	}

	if (extended) {
		mu_names = sprintf("mu%d", 1:L)
		return(c(mu_names, w_names, snames))
	} else {
		return(c(w_names, snames))
	}
}

#' @name cmm_reg
#' @export
print.cmm_reg = function(x, ...)
{
	coefs_ext = par2vec(x$par)
	idx = setdiff(seq_along(coefs_ext), 1:x$L)
	coefs = coefs_ext[idx]

	V = vcov(x)
	se = sqrt(diag(V))
	zval = coefs / se
	pval = 2 * pnorm(abs(zval), lower.tail = FALSE)

	dat = data.frame(
		estimate = coefs,
		se = se,
		zval = zval,
		pval = pval
	)
	rownames(dat) = get_par_names(x)

	printf("Fit for CMM model:\n")
	printf("--- Parameter Estimates ---\n")
	print(dat)

	printf("--\n")
	printf("Iterations: %d   ", x$iter)
	printf("Elapsed Sec: %0.4f   ", x$elapsed_sec)
	printf("Tolerance: %g\n", x$tol)
	printf("Baseline category: %d   ", x$base)
	printf("LogLik: %0.4f   ", logLik(x))
	printf("AIC: %0.4f\n", AIC(x))
}
