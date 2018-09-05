
#' convert a VCF into an 012,-1 matrix and read_depth bin matrix for estimation
#'
#' @param v a vcfR object into which a VCF file has been read
#' @param DF Field to use for obtaining total read depth.  Currently
#' all that is available is DP.
#' @param minBin minimum number of observations for each read depth bin
#' @return This sends back a list with \code{mat012}: the 012,-1 matrix of genotypes, and
#' \code{dp_bins_list}: the list returned by bin_depths().
#' @export
#' @keywords internal
#' @examples
#' pv <- prep_vcf_for_est_m_rd(lobster_buz_2000, minBin = 1000)
prep_vcf_for_est_m_rd <- function(v, DF = "DP", minBin) {

  # check to make sure the field named in DF exists.
  #vt <- vcfR2tidy(v, info_only = TRUE)
  #if(!(DF %in% vt$meta$ID)) {
  #  stop("The tag ", DF, " does not appear to be in the data...")
  #}

  # extract the matrix as an 012 file
  dgt <- vcfR::extract.gt(v, element = "GT")

  d012 <- make_it_012(dgt)
  dimnames(d012) <- dimnames(dgt)

  # now get the read depth matrix
  dp <- vcfR::extract.gt(v, element = "DP")
  storage.mode(dp) <- "integer"
  dp[d012 == -1] <- NA

  # OK, now dp is a matrix of read depths with NAs where the genotype was unobserved

  # now we need to bin those depths up
  bd <- bin_depths(D = dp, S = minBin)
  bd$dp_bins <- t(bd$dp_bins)


  list(mat012 = t(d012),
       dp_bins_list = bd)
}

