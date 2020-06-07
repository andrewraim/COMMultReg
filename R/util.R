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
