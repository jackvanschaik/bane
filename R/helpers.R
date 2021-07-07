#' Sigmoid
#'
#' @param x a real number
#'
#' @return `1/(1 + exp(-x))`
sig <- function(x) {
    1/(1 + exp(-x))
}

#' Inverse Sigmoid
#'
#' Inverse of `1/(1 + exp(-x))`
#'
#' @param x a real number between 0 and 1
#'
#' @return `log(x/(1-x))`
invsig <- function(x) {
    log(x/(1-x))
}
