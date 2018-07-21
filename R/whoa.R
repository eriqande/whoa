#' whoa: Evaluation of genotyping error in genotype by sequencing data
#'
#' The name is found in the capitals here:  \code{W}here's my \code{H}eterozygotes?! \code{O}bservations
#' of genotyping \code{A}ccuracy.
#'
#'
#' @section the \code{whoa} main high-level functions:
#'
#' Fill in.
#'
#'
#' @section genetic data formats:
#'
#' Fill in.
#'
#' @section example data:
#'
#' Fill in.
#'
#' @docType package
#' @name whoa
#' @importFrom Rcpp evalCpp
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr summarise
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#' @importFrom stats setNames
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @useDynLib whoa
NULL


# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if(getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "Freq",
      "dp",
      "bin"
    )
  )
}


