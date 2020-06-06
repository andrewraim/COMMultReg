#' Conway Maxwell Multinomial Regression
#' 
#' Functions for CMM Regression.
#' 
#' @param y 
#' @param m 
#' @param X 
#' @param W 
#' @param base 
#' @param par_init
#' @param control 
#' 
#' @return
#'
#' @examples
#' @name cmm
NULL

#' @export
cmm_reg = function(y, m, X, W, base = 1, par_init = NULL, control = NULL)
{
	dat_xform = transform_data(y, m, X, W)
	k = ncol(y)
	p_x = ncol(X)
	p_w = ncol(W)
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
		base = base, tol = control$tol,
		xnames = colnames(X), wnames = colnames(W))
	out$loglik_poisson = out$loglik
	out$base = base
	out$loglik = loglik_cmm_xform(out$par, dat_xform, base = base)
	return(out)
}

#' @export
cmm_reg_control = function(tol = 1e-8, verbose = FALSE, max_iter = 200)
{
	list(tol = tol, verbose = verbose, max_iter = max_iter)
}

#' @export
par2vec = function(par)
{
	as.numeric(c(par$mu, par$gamma, as.numeric(par$beta)))
}

#' @export
vec2par = function(vartheta, L, p_x, p_w, k)
{
	list(
		mu = vartheta[1:L],
		gamma = vartheta[L + 1:p_w],
		beta = matrix(vartheta[L + p_w + seq_len(p_x*(k-1))], p_x, k-1)
	)
}

#' @export
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

#' @export
loglik_cmm_xform = function(par, dat_xform, base = 1)
{
	L = length(dat_xform)
	vartheta = par2vec(par)
	qq = length(vartheta)
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
				svartheta = sum(s * vartheta)
				ll = ll + dat$freqs[i] * svartheta
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

#' @export
newton_raphson = function(par_init, dat_xform, base = 1, xnames = NULL,
	wnames = NULL, tol = 1e-6, max_iter = 100, verbose = FALSE)
{
	L = length(dat_xform)
	p_x = nrow(par_init$beta)
	p_w = length(par_init$gamma)
	k = ncol(par_init$beta) + 1

	delta = Inf
	vartheta = par2vec(par_init)
	iter = 0
	st = Sys.time()

	if (verbose) {
		out = loglik_score_fim_cmm(par_init, dat_xform, base = base)
		logger("iter 0: Poisson loglik: %g\n", out$loglik)
		print(vartheta[-(1:L)])
	}

	while (delta > tol && iter < max_iter) {
		iter = iter + 1
		vartheta_old = vartheta

		par = vec2par(vartheta, L, p_x, p_w, k)
		out = loglik_score_fim_cmm(par, dat_xform, base = base)
		move = tryCatch({
			solve(out$fim, out$score)
		}, error = function(e) {
			warning("FIM was close to singular, trying pseudo-inverse instead...")
			pinv(out$fim) %*% out$score
		})
		vartheta = vartheta + move
		delta = sum(abs(vartheta - vartheta_old))

		if (verbose) {
			logger("iter %d: Poisson loglik: %g delta: %g\n", iter, out$loglik, delta)
			print(vartheta[-(1:L)])
		}
	}

	elapsed_sec = as.numeric(Sys.time() - st, units = "secs")

	par = vec2par(vartheta, L, p_x, p_w, k)
	colnames(par$beta) = sprintf("cat%d", (1:k)[-base])
	if (is.null(xnames)) {
		rownames(par$beta) = sprintf("X%d", 1:p_x)
	} else if (length(xnames) != p_x) {
		warning("length of xnames does not match ncol(X)")
	} else {
		rownames(par$beta) = xnames
	}

	ret = list(par = par, loglik = out$loglik, score = out$score,
		fim = out$fim, tol = delta, iter = iter, L = L,
		xnames = xnames, wnames = wnames, elapsed_sec = elapsed_sec)
	class(ret) = "cmm_reg"
	return(ret)
}

#' @export
logLik.cmm_reg = function(object, ...)
{
	object$loglik
}

#' @export
AIC.cmm_reg = function(object, ..., k = 2)
{
	p = length(object$par$gamma) + length(object$par$beta)
	-2*logLik(object) + k*p
}

#' @export
coef.cmm_reg = function(object, ...)
{
	vartheta = par2vec(object$par)
	idx = setdiff(1:length(vartheta), 1:object$L)
	vartheta[idx]
}

#' @export
vcov.cmm_reg = function(object, extended = FALSE, ...)
{
	vartheta = par2vec(object$par)	
	idx = setdiff(1:length(vartheta), 1:object$L)
	V = tryCatch({
		solve(object$fim)
	}, error = function(e) {
		warning("FIM was close to singular, trying pseudo-inverse instead...")
		pinv(object$fim)
	})
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

#' @export
print.cmm_reg = function(x, ...)
{
	vartheta = par2vec(x$par)
	idx = setdiff(1:length(vartheta), 1:x$L)
	V = vcov(x, extended = TRUE)
	se = sqrt(diag(V))
	zval = vartheta / se
	pval = 2 * (1 - pnorm(abs(zval)))

	dat = data.frame(
		estimate = vartheta[idx],
		se = se[idx],
		zval = zval[idx],
		pval = pval[idx]
	)

	k = ncol(x$par$beta) + 1
	p_x = nrow(x$par$beta)
	p_w = nrow(x$par$gamma)

	xnames = x$xnames
	cats = setdiff(1:k, x$base)
	if (is.null(xnames)) {
		label_var = rep(1:p_x, times = k-1)
		label_dim = rep(cats, each = p_x)
		snames = sprintf("%d:beta%d", label_dim, label_var)
	} else {
		stopifnot(length(xnames) == nrow(x$par$beta))
		label_var = rep(xnames, times = k-1)
		label_dim = rep(cats, each = p_x)
		snames = sprintf("%d:%s", label_dim, xnames)
	}

	wnames = x$wnames
	if (is.null(wnames)) {
		wnames = sprintf("gamma%d", 1:length(x$par$gamma))
	} else {
		stopifnot(length(wnames) == length(x$par$gamma))
	}

	rownames(dat) = c(wnames, snames)

	printf("Fit for CMM model:\n")
	printf("--- Parameter Estimates ---\n")
	print(dat)

	printf("--\n")
	printf("Iterations: %d   ", x$iter)
	printf("Elapsed Sec: %0.4f   ", x$elapsed_sec)
	printf("Tolerance: %g\n", x$tol)
	printf("LogLik: %0.4f   ", logLik(x))
	printf("AIC: %0.4f\n", AIC(x))
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
