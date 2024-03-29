#' Conway-Maxwell-Multinomial Regression
#' 
#' @details
#' This package provides basic drawing and density functions for the
#' Conway-Maxwell-Multinomial (CMM) distribution, as well as functions to
#' carry out maximum likelihood estimation (MLE) for a regression model
#' based on CMM. The MLE method was inspired by Altham and Hankin (2012).
#' Further details can be found in Morris, Raim, and Sellers (2020).
#' 
#' Computationally intensive sections are carried out in C++ with the help
#' of Rcpp and RcppArmadillo. Logic to iterate through the multinomial sample
#' space, without pre-generating the elements, is from the function
#' \code{gsl_multiset_next} in the GNU GSL library (Galassi et al).
#' 
#' Note that some of the functions in this package currently iterate through
#' the multinomial sample space based on \eqn{m} trials and \eqn{k} categories.
#' These functions are suitable for situations in which \eqn{m+k-1 \choose m}
#  is not very large for any particular observation.
#' 
#' @useDynLib COMMultReg, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats AIC coef logLik model.frame model.matrix model.offset 
#' model.response model.weights pnorm vcov
#' @importFrom combinat xsimplex
#' 
#' @name COMMultReg-package
#' 
#' @references 
#' Pat M. E. Altham and Robin K. S. Hankin (2012). Multivariate generalizations
#' of the multiplicative binomial distribution: Introducing the MM package.
#' Journal of Statistical Software 46.
#' 
#' Darcy Steeg Morris, Andrew M. Raim, and Kimberly F. Sellers (2020).
#' A Conway-Maxwell-multinomial distribution for Flexible modeling
#' of clustered categorical data. Journal of Multivariate Analysis, 179:104651.
#' Preprint: https://arxiv.org/abs/1911.02131.
#' 
#' Dirk Eddelbuettel (2013) Seamless R and C++ Integration with Rcpp.
#' Springer, New York. ISBN
#' 
#' Dirk Eddelbuettel, Conrad Sanderson (2014). RcppArmadillo: Accelerating
#' R with high-performance C++ linear algebra. Computational Statistics and
#' Data Analysis, Volume 71, March 2014, pages 1054-1063.
#' 
#' M. Galassi et al, GNU Scientific Library Reference Manual, 3rd Edition,
#' ISBN 0954612078.
NULL
