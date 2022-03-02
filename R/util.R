#' printf
#' 
#' Print a message with similar arguments as \code{sprintf}.
#' 
#' @param msg format string.
#' @param ... values to be passed into msg.
#'
#' @export
printf = function(msg, ...)
{
	cat(sprintf(msg, ...))
}

#' logger
#' 
#' Print a message proceeded by a timestamp, with similar arguments as \code{sprintf}.
#' 
#' @param msg format string.
#' @param ... values to be passed into msg.
#'
#' @export
logger = function(msg, ...)
{
	sys.time = as.character(Sys.time())
	cat(sys.time, "-", sprintf(msg, ...))
}

#' logger
#' 
#' Inverse multinomial logit transformation.
#' 
#' @param x a \eqn{p-1} dimensional vector.
#' @param base Index of the baseline category.
#' 
#' @returns a vector in the \eqn{p} dimensional probability simplex.
#'
#' @export
inv_mlogit = function(x, base = 1) {
	p = numeric(length(x)+1)
	z = exp(x)
	p[base] = 1 / (1 + sum(z))
	p[-base] = z * 1 / (1 + sum(z))
	return(p)
}

#' logger
#' 
#' Inverse multinomial logit transformation.
#' 
#' @param p a vector in the \eqn{p} dimensional probability simplex.
#' @param base Index of the baseline category.
#' 
#' @returns a \eqn{p-1} dimensional vector.
#'
#' @export
mlogit = function(p, base = 1) {
	log(p[-base] / p[base])
}
