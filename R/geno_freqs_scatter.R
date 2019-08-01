

#' return a `ggplot2' plot object of observed and expected genotype freqs
#'
#' @param gfc a tibble like that created by exp_and_obs_geno_freqs()
#' @param alpha the transparency (alpha) parameter to apply to the points
#' in the scatterplot. Default is 0.2.
#' @param max_plot_loci By default this plots only 500 loci, sampled
#' randomly, to keep `ggplot2' taking forever to plot, for example, 100K
#' points.  If you want to plot all the points, set this to a number
#' larger than the number of single nucleotide polymorphisms (SNPs) in the data set.
#' @export
#' @examples
#' # get the expected and observed geno freqs
#' gfreqs <- exp_and_obs_geno_freqs(lobster_buz_2000)
#' g <- geno_freqs_scatter(gfreqs)
#'
#' # now g is a 'ggplot2' object.
geno_freqs_scatter <- function(gfc, alpha = 0.2, max_plot_loci = 500) {

  snps <- unique(gfc$snp)
  if(length(snps) > max_plot_loci) {
    ss <- sample(x = snps, size = max_plot_loci, replace = FALSE)
  } else {
    ss <- snps
  }
  g <- ggplot2::ggplot(gfc %>% dplyr::filter(snp %in% ss) , ggplot2::aes(x = p_exp, y = p_obs, colour = geno)) +
    ggplot2::geom_jitter(alpha = alpha, position = ggplot2::position_jitter(width = 0.01, height = 0.01)) +
    ggplot2::facet_wrap(~ geno, nrow = 1) +
    ggplot2::geom_polygon(data = geno_freq_boundaries(), fill = NA, linetype = "dashed", colour = "black") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "solid")

  g

}
