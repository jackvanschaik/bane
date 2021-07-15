#' Estimate Normal Distribution of Inverse Sigmoid Transformed Beta
#'
#' Given beta shape parameters A and B, this will estimate the mean
#' and variance of the inverser sigmoid transformed distribution via
#' a normal distribution.
#'
#' This is useful for choosing prior values. Prior values should be given in
#' terms of log odds, but sometimes proportions are more intuitively appealing.
#' This function will estimate the log odds distibution given the corresponding
#' beta distribution of a proprtion.
#'
#'
#' @param A The alpha shape parameter
#' @param B The beta shape parameter
#' @param size The sample size
#'
#' @return a list of values of interest
#' @export
estimate_normal <- function(A, B, size = 100000) {
    A <- as.numeric(A)
    B <- as.numeric(B)
    x_b <- stats::rbeta(size, A, B)
    x_n <- invsig(x_b)
    list(
        sample = data.frame(beta = x_b, norm = x_n),
        mean_b = mean(x_b),
        var_b = stats::var(x_b),
        mean_n = mean(x_n),
        var_n = stats::var(x_n)
    )
}
