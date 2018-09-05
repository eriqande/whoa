

#' get posterior estimates for m from different read depth categories
#'
#' This function calls internal C++ routines that perform Markov
#' chain Monte Carlo to sample from the posterior distrubtion of the
#' heterozygote miscall rate for each read depth category.
#'
#' The read depth bins are
#' determined by giving the function minBin---the minimum number of observations desired for each read
#' depth bin.  The function then breaks the observations into bins so that each read depth bin
#' has at least minBin observations.
#'
#' Note that if you want to estimate the heterozygote miscall rate \bold{overall} (i.e., not conditioning each
#' estimate on a read depth bin), then simply give a very large number (larger than the number of
#' markers times the number of individuals) for minBin.  For example, you could use a number like 1e15
#' for minBin. As a consequence, all the genotypes will be put into a single read depth bin.
#' @param v a vcfR object holding the information from a VCF file with the genotype and depth data.
#' If you are going to be it down by read depth categories, the VCF file must have a DP field for
#' every genotype.
#' @param indivs a character vector holding
#' the names of the individuals from the VCF file to include in the analysis.  By default
#' this is NULL, in which case everyone from the file is included.
#' @param minBin minimum number of observations for each read depth bin.  If you have 10K markers
#' and 50 individuals then you have about 500,000 genotypes to play with.  Requiring bins with at
#' least 5,000 genotypes in them will give you less than 100 bins.  You can play around with this
#' number to get the right number of bins.  The algorithm breaks the read depths up into bins that
#' have at least minBin genotypes in them.
#' @param init_m the initial value of the heterozygote miscall rate assumed for each read depth bin.
#' By default this is 0.1.
#' @param num_sweeps the number of sweeps of the MCMC algorithm to do.  Default is 500, which is a
#' little on the short side.  Run multiple times and make sure the values obtained are similar across
#' runs to assess convergence.
#' @param burn_in how many sweeps from the beginning of the chain to discard when computing the
#' posterior means and quantiles.  Default is 100.  Note that full traces of the visited m values
#' for every read depth bin are returned as well, so that the behavior of the chain in those early
#' steps can be investigated.
#' @return A list with six components:
#' \describe{
#'   \item{m_posteriors}{A tibble with 6 columns: bin = the index of the read depth bin;
#'   mean = the posterior mean estimate of the the heterozygote miscall rate in that bin;
#'   lo95 = the low endpoint of the 95% credible interval; hi95 = the high endpoint of the
#'   95% credible interval; total_n = the total number of genotypes in the read depth bin;
#'   and mean_dp = the mean read depth in the bin. }
#'   \item{m_traces}{A tibble with all the values visited for m for every read depth bin.  This tibble
#'   has three columns: bin = the index of the read depth bin; sweep = the sweep number, value = the value
#'   of m for that read depth bin in that particular sweep. }
#'   \item{dp_summary}{A tibble summarizing how many genotypes of different read depths
#'   appear in each bin.}
#'   \item{bin_stats}{A tibble with a different summary of the read-depth bins.}
#'   \item{num_sweeps}{Number of MCMC sweeps used.}
#'   \item{burn_in}{Number of sweeps discarded as burn in.}
#' }
#' @export
#' @examples
#' # Shorter run than recommended (for quick example...)
#' im <- infer_m(lobster_buz_2000, minBin = 1000, num_sweeps = 100, burn_in = 20)
infer_m <- function(v, minBin, indivs = NULL, init_m = 0.1, num_sweeps = 500, burn_in = 100) {

  message("Preparing data structures for MCMC")

  # deal with the indivs
  if(is.null(indivs)) {
    indivs <- colnames(v@gt)[-1]
  }

  # check to make sure that all of indivs are included in the VCF
  wrongos <- setdiff(indivs, colnames(v@gt)[-1])
  if(length(wrongos) > 0) {
    stop("Requesting indivs that are not in v: ", paste(wrongos, collapse = ", "))
  }

  # extract just those indivs.  This assumes that vcfR has named the first column FORMAT
  v@gt <- v@gt[,c("FORMAT", indivs)]




  # prep it
  p <- prep_vcf_for_est_m_rd(v, "DP", minBin)


  message("Running MCMC")
  # estimate it
  b <- estimate_m_rd(Y = p$mat012,
                     R = p$dp_bins_list$dp_bins,
                     init_m = init_m,
                     num_cats = p$dp_bins_list$num_cats,
                     p_prior = c(0.5, 0.5),
                     m_prior = c(0.5, 0.5),
                     num_reps = num_sweeps)

  message("Tidying output")
  # tidy the output
  res <- tidy_m_ests(E = b,
                     S = p$dp_bins_list$bin_stats,
                     TB = p$dp_bins_list$tidy_bins,
                     burn = burn_in)

  res
}
