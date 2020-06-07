#' Conway Maxwell Multinomial Regression
#' 
#' Functions for CMM Regression.
#' 
#' @param formula_x
#' @param formula_w
#' @param y TBD
#' @param X TBD
#' @param W TBD
#' @param base TBD
#' @param extended TBD
#' @param par_init TBD
#' @param control TBD
#' 
#' @return
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
#' changed to category \code{a} by specifying the \code{base = a} argument.
#'
#' @examples
#' print("TBD")
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

	# TBD: Do something better with par_init
	cmm_reg_raw(y, X, W, par_init = NULL, control = control)
}

#' @name cmm_reg
#' @export
cmm_reg_raw = function(y, X, W, par_init = NULL, control = NULL)
{
	m = rowSums(y)
	k = ncol(y)
	p_x = ncol(X)
	p_w = ncol(W)

	dat_xform = transform_data(y, m, X, W)
	L = length(dat_xform)

	if (is.null(par_init)) {
		par_init = list(
			mu = numeric(L),
			beta = matrix(0, p_x, k-1),
			gamma = numeric(p_w)
		)
	} else {
		check_par(par_init, dat_xform)
	}

	if (is.null(control)) {
		control = cmm_reg_control()
	}

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
	out$xnames = colnames(X)
	colnames(out$par$beta) = sprintf("cat%d", (1:k)[-control$base])
	if (is.null(out$xnames)) {
		rownames(out$par$beta) = sprintf("X%d", 1:p_x)
	} else if (length(out$xnames) != p_x) {
		warning("length of xnames does not match ncol(X)")
	} else {
		rownames(out$par$beta) = out$xnames
	}

	out$wnames = colnames(W)
	names(out$par$gamma) = out$wnames

	return(out)
}

#' @name cmm_reg
#' @export
cmm_reg_control = function(base = 1, tol = 1e-8, verbose = FALSE, max_iter = 200)
{
	list(base = base, tol = tol, verbose = verbose, max_iter = max_iter)
}

check_par = function(par, dat_xform)
{
	stopifnot(!is.null(par_init[["mu"]]))
	stopifnot(!is.null(par_init[["beta"]]))
	stopifnot(!is.null(par_init[["gamma"]]))
	mu = par_init[["mu"]]
	beta = par_init[["beta"]]
	gamma = par_init[["gamma"]]

	L = length(dat_xform)
	stopifnot(L == length(mu))
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

transform_data = function(y, m, X, W, tol = 1e-8)
{
	n = nrow(y)
	k = ncol(y)
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

loglik_cmm_xform = function(par, dat_xform, base = 1)
{
	L = length(dat_xform)
	coefs = par2vec(par)
	qq = length(coefs)
	ll = 0

	if (TRUE) {
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
	} else {
		# This part should no longer be needed. Remove it?
		# At least we can use the collapsed form of the data to avoid extra walks
		# through the multinomial sample space
		for (l in 1:L) {
			dat = dat_xform[[l]]
			for (i in 1:nrow(dat$z_obs)) {
				ll = ll + dat$freqs[i] * d_cmm(
					x = dat$z_obs[i,],
					m = dat$m,
					p = t(dat$x) %*% inv_mlogit(par$beta),
					nu = t(dat$w) %*% par$gamma,
					take_log = TRUE, normalize = FALSE)
			}
			ll = ll - sum(dat$freqs) * d_cmm_normconst(m = dat$m,
					p = t(dat$x) %*% inv_mlogit(par$beta),
					nu = t(dat$w) %*% par$gamma,
					take_log = TRUE)
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
		out = loglik_score_fim_cmm(par_init, dat_xform, base = base)
		logger("iter 0: Poisson loglik: %g\n", out$loglik)
		print(coefs[-(1:L)])
	}

	while (delta > tol && iter < max_iter) {
		iter = iter + 1
		coefs_old = coefs

		par = vec2par(coefs, L, p_x, p_w, k)
		out = loglik_score_fim_cmm(par, dat_xform, base = base)
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

#' @name cmm_reg
#' @export
get_par_names = function(object, extended = FALSE, ...)
{
	par = coef(object, extended = TRUE)
	L = length(par$mu)
	k = ncol(par$beta) + 1
	p_x = nrow(par$beta)
	p_w = length(par$gamma)

	xnames = object$xnames
	cats = setdiff(1:k, object$base)
	if (is.null(xnames)) {
		label_var = rep(1:p_x, times = k-1)
		label_dim = rep(cats, each = p_x)
		snames = sprintf("%d:beta%d", label_dim, label_var)
	} else {
		stopifnot(length(xnames) == nrow(par$beta))
		label_var = rep(xnames, times = k-1)
		label_dim = rep(cats, each = p_x)
		snames = sprintf("%d:%s", label_dim, xnames)
	}

	wnames = object$wnames
	if (is.null(wnames)) {
		wnames = sprintf("gamma%d", 1:length(par$gamma))
	} else {
		stopifnot(length(wnames) == length(par$gamma))
	}

	if (extended) {
		munames = sprintf("mu%d", 1:L)
		return(c(munames, wnames, snames))
	} else {
		return(c(wnames, snames))
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
	pval = 2 * (1 - pnorm(abs(zval)))

	dat = data.frame(
		estimate = coefs[idx],
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
