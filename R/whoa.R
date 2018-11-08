#' whoa: Evaluation of genotyping error in genotype by sequencing data
#'
#' The name is found in the capitals here:  \code{W}here's my \code{H}eterozygotes?! \code{O}bservations
#' of genotyping \code{A}ccuracy.
#'
#'
#' @section the \code{whoa} main high-level functions:
#'
#' The main function in the package whoa is \code{\link{infer_m}}.
#' This function infers the heterozygote miscall rate (the rate at
#' which true heterozygotes have been miscalled as homozygotes in
#' genotype-by-sequencing data) for calls made upon genotypes falling
#' within different read depth categories.
#'
#' The output from  \code{\link{infer_m}} is easily plotted by passing the
#' m_posteriors component of the output list from \code{infer_m} into
#' \code{\link{posteriors_plot}}.
#'
#'
#' @section example data:
#'
#' The package comes with a small data set, \code{\link{lobster_buz_2000}}, which
#' was read in from a VCF file and is now stored in the package as a vcfR object.
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
      "bin",
      "mean_dp",
      "lo95",
      "hi95",
      "total_n",
      "0",
      "1",
      "2",
      "geno",
      "n_exp",
      "n_obs",
      "ntot",
      "p_exp",
      "p_obs",
      "snp",
      "z_score",
      "INDIVIDUALS",
      "ID"
    )
  )
}


