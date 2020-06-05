# We needed some workarounds to get the package to work. See:
# - https://stackoverflow.com/questions/36605955/c-function-not-available
# - https://cran.r-project.org/web/packages/roxygen2/vignettes/namespace.html

#' @useDynLib COMMultReg, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL
