#' BaNE examples
#'
#' Run one of a few examples using a built-in dataset
#'
#' @param example One of "all_conditions" or "two_conditions". View source code for model specifications.
#'
#' @return A plot of posterior subcohort proportions
#' @export
bane_example <- function(example = "two_conditions") {
    if (example == "all_conditions") {
        ex <- Bane$new(
            indep = c("old", "male"),
            dep = list(
                white = c("old", "male"),
                hyper = c("white", "old")
            ),
            mu = c(0.4, 0.49),
            lm = c(0.8, 0.25),
            data = bane::synthetic_ehr
        )
        ex$run_chain(Iterations = 250, Status = 25, Thinning = 5, L = 1, eps = 0.02)
        return(plot(ex$post_subs()$ggplot))
    }
    else if (example == "two_conditions") {
        ex <- Bane$new(
            indep = c("male"),
            dep = list(
                white = c("male")
            ),
            mu = c(0.5),
            lm = c(0.8),
            data = bane::synthetic_ehr
        )
        ex$run_chain(Iterations = 250, Status = 25, Thinning = 5, L = 1, eps = 0.02)
        return(plot(ex$post_subs()$ggplot))
    }
    else {
        stop("Try example = 'two_conditions' or 'all_conditions'")
    }
}
