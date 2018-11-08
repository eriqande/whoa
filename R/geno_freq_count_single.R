#' from on 012 file compute expected (assuming HW equilbrium) and observed genotype counts
#'
#' This is an internal function.
#' @param g012 an 012 matrix, with indivs in rows and loci in columns.
#' Missing data can be -1 or NA. The matrix must
#' have colnames which are the locus names.
#' @export
#' @keywords internal
#' @examples
#' # get an 012 matrix from the lobster data
#' tmp <- t(vcfR::extract.gt(lobster_buz_2000, element = "GT"))
#' locnames <- colnames(tmp)
#' g <- make_it_012(tmp)
#' colnames(g) <- locnames # put these back on since make_it_012 removes them
#' gf <- geno_freq_calc_single(g)
#' @return See information for return value of \code{\link{exp_and_obs_geno_freqs}}.
#' @importFrom tidyr gather

geno_freq_calc_single <- function(g012) {
  g012[g012 == -1] <- NA

  gc <- 2 * colSums(!is.na(g012), na.rm = TRUE)  # number of non-missing gene copies
  p <- colSums(g012, na.rm = TRUE) / gc
  q <- 1.0 - p

  n0 <- colSums(g012 == 0, na.rm = TRUE)
  n1 <- colSums(g012 == 1, na.rm = TRUE)
  n2 <- colSums(g012 == 2, na.rm = TRUE)
  n <- n0 + n1 + n2

  tibble::tibble(snp = colnames(g012),
                 p = p,
                 ntot = n,
                 `0` = n0,
                 `1` = n1,
                 `2` = n2
  ) %>%
    tidyr::gather(key = "geno", value = "n_obs", `0`, `1`, `2`) %>%
    dplyr::arrange(snp, geno) %>%
    dplyr::mutate(p_exp = ifelse(geno == 0, (1 - p)^2,
                                 ifelse(geno == 1, 2 * p * (1 - p),
                                        ifelse(geno == 2, p^2, NA))),
                  n_exp = ntot * p_exp,
                  p_obs = n_obs / ntot,
                  z_score = (n_obs - n_exp) / sqrt(ntot * p_exp * (1 - p_exp))) %>%
    dplyr::select(snp, p, ntot, geno, p_exp, p_obs, n_exp, n_obs, z_score)


}
