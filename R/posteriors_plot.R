

#' @title Plot the posterior estimates for heterozyote miscall rates
#'
#' @description This just returns a ggplot object that plots the read depth bins on the x-axis
#' and the posterior mean m estimates (and credible intervals) on the y-axis, and
#' depicts the number of genotypes in each read depth bin using color.

#' @param P the tibble that is the m_posteriors component of \code{\link{infer_m}}
#' @return a ggplot2 object.
#' @export
#' @examples
#' # get something to plot (short run for example)
#' im <- infer_m(lobster_buz_2000, minBin = 1000, num_sweeps = 100, burn_in = 20)
#'
#' # then plot it
#' g <- posteriors_plot(im$m_posteriors)
#'
#' # now g is a ggplot object
posteriors_plot <- function(P) {
  ggplot2::ggplot(P) +
    ggplot2::geom_line(ggplot2::aes(x = mean_dp, y = mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(x = mean_dp, ymin = lo95, ymax = hi95), fill = "pink", alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = mean_dp, y = mean, colour = total_n)) +
    viridis::scale_colour_viridis() +
    ggplot2::xlab("Mean read depth of bin") +
    ggplot2::ylab("Posterior mean estimate of miscall rate")
}
