

#' bin read depths of SNPs into categories having at least S observations
#'
#' @param D a matrix of read depths.  Rows are individuals, columns are SNPs.  Cells where data are missing
#' in the genotype matrix must be denoted as NA
#' @param S the min number of observations to have in each bin
#' @return This returns a list with two components.  \code{dp_bins} is a matrix of the same
#' shape as D with the bin categories (as 1, 2, ...) and -1 for this cells
#' corresponding to missing genotypes.  \code{num_cats} is the number of depth bins.
#' \code{tidy_bins} is a long format description of the bins.
#' \code{bin_stats} is a tibble giving summary information about the read depth bins which
#' is useful for plotting things, etc.
#' @export
#' @keywords internal
#' @examples
#'
#' # get a matrix of read depths and make it an integer matrix
#' depths <- vcfR::extract.gt(lobster_buz_2000, element = "DP")
#' storage.mode(depths) <- "integer"
#'
#' # get a character matrix of genotypes, so we can figure out which
#' # are missing and mask those from depths
#' genos <- vcfR::extract.gt(lobster_buz_2000, element = "GT")
#'
#' # make missing in depths if missing in genos
#' depths[is.na(genos)] <- NA
#'
#' # bin the read depths into bins with at least 1000 observations in each bin
#' bins <- bin_depths(depths, 1000)
bin_depths <- function(D, S) {

  # first, count things up and get them into a reasonable format
  dcnts <- tibble::as_tibble(as.data.frame(table(D))) %>%
    dplyr::mutate(D = as.integer(as.character(D))) %>%
    setNames(c("D", "Freq")) %>%
    dplyr::rename(dp = D,
           n = Freq) %>%
    dplyr::arrange(dp) # just to make sure it is in increasing order

  # now we lump them up.  We just use a for loop for this
  idx <- 1
  sum <- 0
  n <- dcnts$n
  res <- rep(NA, length(n))
  for(i in seq_along(n)) {
    res[i] <- idx
    sum <- sum + n[i]
    if(sum > S) {
      idx <- idx + 1
      sum <- 0
    }
  }

  # at the end of that, we take the very last bin and just merge is with the
  # preceding one, so that it doesn't have just a small number in it.
  # but we only do this if there is more than one read depth bin!
  if(max(res) > 1) {
    res[res == max(res)] <- max(res) - 1
  }

  # and we add it on there with a mutate which will bark an error is something has gone awry
  tidy_bins <- dcnts %>%
    dplyr::mutate(bin = as.integer(res))

  # now we have to cut the original matrix up into these categories
  left_endpoints <- tidy_bins %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(ends = max(dp) + 0.2)

  cut_vec <- c(0, left_endpoints$ends)


  dp_bins <- D
  dp_bins[] <- cut(D, breaks = cut_vec)

  dp_bins[is.na(dp_bins)] <- -1L

  # finally summarize the bin stats
  bin_stats <- tidy_bins %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(total_n = sum(n),
              mean_dp = sum(dp * n) / sum(n))

  list(dp_bins = dp_bins,
       num_cats = max(tidy_bins$bin),
       tidy_bins = tidy_bins,
       bin_stats = bin_stats
  )
}
