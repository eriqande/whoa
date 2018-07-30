
#' tidy up the estimate_m_rd output into something you can plot
#' @param E the list returned by the estimation function
#' @param S the bin stats
#' @param TB the tidy bins
#' @param burn how much for burn in?
#' @keywords internal
#' @export
tidy_m_ests <- function(E, S, TB, burn = 50) {

  # get the posteriors and the quantiles
  posts <- tibble::tibble(bin = 1:nrow(E$Mtrace),
         mean = rowMeans(E$Mtrace[, -(1:burn), drop = FALSE]),
         lo95 = apply(E$Mtrace[, -(1:burn), drop = FALSE], 1, quantile, probs = 0.05),
         hi95 = apply(E$Mtrace[, -(1:burn), drop = FALSE], 1, quantile, probs = 0.95)
  ) %>%
    left_join(S, by = "bin")

  # now, also get the traces for m values
  trace <- tibble::tibble(
    bin = rep(1:nrow(E$Mtrace), length.out = length(E$Mtrace)),
    sweep = rep(1:ncol(E$Mtrace), each = nrow(E$Mtrace)),
    value = as.numeric(E$Mtrace)
  )

  list(m_posteriors = posts,
       m_traces = trace,
       dp_summary = TB,
       bin_stats = S,
       num_sweeps = ncol(E$Mtrace),
       burn_in = burn
       )
}

