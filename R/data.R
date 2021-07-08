#' Synthetic EHR Test Data
#'
#' A synthetic dataset containing four condition flag, meant to test BaNE.
#' There is one row per synthetic patient. Original data was generated via
#' Synthea, then transformed to flags.
#'
#' @format A data frame with 5995 rows and 4 variables
#' \describe{
#' \item{old}{1 if patient was born before 1970, 0 otherwise}
#' \item{male}{1 if patient's sex is male, 0 otherwise }
#' \item{white}{1 if patient's race is white, 0 otherwise}
#' \item{hyper}{1 if patient had a hypertension diagnosis on record, 0 otherwise}
#' }
"synthetic_ehr"
