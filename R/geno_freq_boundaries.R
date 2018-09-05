
#' function to return a tibble with the min/max values of genotype freqs possible
#'
#' These mins and maxes occur because the genotypes are used to estimate
#' the allele frequencies.  The output of this function is useful for
#' putting boundaries on plots.
#' @return a tibble
#' @export
#' @keywords internal
#' @examples
#' gfb <- geno_freq_boundaries()
geno_freq_boundaries <- function() {
  # first do it for the homozygote category
  Poe <- seq(0,1, by = 0.005)
  phat <- 1 - sqrt(Poe)
  minPo <- pmax(0, 1 - 2 * phat)
  maxPo <- 1 - phat

  # these are the values for the two homozygote categories
  homo_tib <- tibble::tibble(p_exp = rep(Poe, 4),
                             p_obs = rep(c(minPo, maxPo), 2),
                             geno = as.character(rep(c(0L, 2L), each = length(Poe) * 2)))

  # now, it should be easy to get the max/min values for heterozygotes.
  # They will occur where one of the homozygotes is min or max.
  P1e <- 2 * phat * (1 - phat)
  maxP1 <- 2 * (1 - phat - minPo)
  minP1 <- 2 * (1 - phat - maxPo)

  het_tib <- tibble::tibble(p_exp = rep(P1e, 2),
                            p_obs = c(minP1, maxP1),
                            geno = "1")

  dplyr::bind_rows(homo_tib, het_tib)

}
