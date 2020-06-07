#' @export
printf = function(msg, ...)
{
	cat(sprintf(msg, ...))
}

#' @export
logger = function(msg, ...)
{
	sys.time = as.character(Sys.time())
	cat(sys.time, "-", sprintf(msg, ...))
}

#' @export
inv_mlogit = function(x, base = 1) {
	p = numeric(length(x)+1)
	z = exp(x)
	p[base] = 1 / (1 + sum(z))
	p[-base] = z * 1 / (1 + sum(z))
	return(p)
}

#' @export
mlogit = function(p, base = 1) {
	log(p[-base] / p[base])
}

#' @export
mult2cat = function(Y) {
	stopifnot(all(sort(unique(as.numeric(Y))) == c(0,1)))
	stopifnot(all(rowSums(Y) == 1))

	n = nrow(Y)
	p = ncol(Y)
	z = integer(n)
	for (j in 1:p) {
		idx = which(Y[,j] == 1)
		z[idx] = j
	}
	return(z)
}

#' @export
cat2mult = function(Y, levels) {
	n = length(Y)
	k = length(levels)
	res = matrix(0, n, k)
	for (j in 1:k) {
		idx = which(Y == levels[j])
		res[idx,j] = 1
	}
	
	colnames(res) = levels
	return(res)
}

