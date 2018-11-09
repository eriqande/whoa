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
#' @importFrom dplyr arrange group_by left_join mutate rename summarise select distinct n_distinct ungroup tally filter
#' @importFrom stats quantile setNames
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr keep flatten_chr
#' @importFrom readr read_tsv write_tsv
#' @importFrom stringi stri_join stri_detect_fixed stri_sub stri_replace_all_fixed
#' @importFrom tidyr gather
#' @importFrom parallel detectCores
# For magrittr: see utils.R

#' @useDynLib whoa
NULL


