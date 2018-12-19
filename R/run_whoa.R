# run_whoa
# wrapper function
# read vcf for whoa

#' @name run_whoa
#' @title Wrapper function that runs whoa
#' @description Wrapper function that runs \href{https://github.com/eriqande/whoa}{whoa} on a vcf.
#' When the strata argument is used, the heterozygote miscall rate is calculated for each strata.

#' @inheritParams read_whoa
#' @param strata (optional) A tab-separated file with 2 columns: \code{INDIVIDUALS} and \code{STRATA}.
#' Using this argument will run the function serially on all the strata and overall.
#' Default: \code{strata = NULL}.
#' @inheritParams exp_and_obs_geno_freqs
#' @inheritParams geno_freqs_scatter
#' @inheritParams infer_m
#' @inheritParams posteriors_plot

#' @param ... (optional) Advance mode that allows to pass further arguments
#' for fine-tuning the function (see details).

#' @return When strata argument is used, a list with:
#' \enumerate{
#' \item \code{gds}: the gds object
#' \item \code{het.miscall.rate.analysis}: a list with \code{gfreqs},
#' \code{plot.freqs}, \code{binned}, \code{plot.post} for each strata.
#' \item \code{random.seed}: the random seed used
#' \item \code{overall.rate}: the overall heterozygote miscall rate (across read depth and strata)
#' }


#' @details
#' \strong{Advance mode, using \emph{dots-dots-dots ...}}
#' \enumerate{
#' \item \code{sample.snps}: (optional, integer) Select a random number of SNPs
#' to speed up analysis, e.g. using \code{sample.snps = 1000}.
#' Default: \code{sample.snps = NULL}, will use all SNPs.
#' \item \code{random.seed} (optional, integer) For reproducibility, set an integer
#' that will be used inside function that requires randomness. With default,
#' a random number is generated.
#' Default: \code{random.seed = NULL}.
#' }

#' @export
#' @rdname run_whoa

#' @examples
#' \dontrun{
#' require(httr)
#' require(ggpubr)
#' library(whoa)
#' # get the vcf and strata files:
#' writeBin(httr::content(httr::GET("https://datadryad.org/bitstream/handle/
#' 10255/dryad.108679/10156-586.recode.vcf?sequence=1"), "raw"),"lobster_data_2015.vcf")
#' writeBin(httr::content(httr::GET("https://datadryad.org/bitstream/handle/
#' 10255/dryad.108679/README.txt?sequence=2"), "raw"),"lobster_strata_2015.tsv")

#' # work on the strata:
#' readr::read_tsv(file = "lobster_strata_2015.tsv",
#'                col_names = c("INDIVIDUALS", "STRATA"),
#'                col_types = "cc") %>%
#'  dplyr::mutate(
#'    STRATA = stringi::stri_sub(str = INDIVIDUALS, from = 1, to = 3),
#'    STRATA = stringi::stri_replace_all_fixed(
#'      str = STRATA,
#'      pattern = c("CAP", "EDN", "SID"),
#'      replacement = c("MAG", "MAG", "DIN"),
#'      vectorize_all = FALSE)
#'  ) %>%
#'  dplyr::filter(STRATA != "TON") %>%
#'  dplyr::arrange(STRATA, INDIVIDUALS) %>%
#'  readr::write_tsv(x = ., path = "lobster_strata_2015.tsv")
#'
#'  # run whoa on all strata
#'  het.miscall.rate.pop <- run_whoa(
#'      data = "lobster_data_2015.vcf", strata = "lobster_strata_2015.tsv",
#'      sample.snps = 2000,
#'      max_plot_loci = 2000,
#'      minBin = 2000)
#'
#'
#' # To view the combined ggplot of the posterior estimates for heterozyote miscall rates
#' fig <- purrr::map(.x = het.miscall.rate.pop$het.miscall.rate.analysis, .f = 4) %>%
#'     ggpubr::ggarrange(plotlist = ., legend = "right", common.legend = TRUE,
#'         labels = names(.), label.x = 0.5, label.y = 0.9)
#' fig
#' }


#' @author Eric Anderson \email{eric.anderson@noaa.gov} and
#' Thierry Gosselin \email{thierrygosselin@@icloud.com}


run_whoa <- function(
  data,
  indivs = NULL,
  strata = NULL,
  prop_indv_required = 0.5,
  prop_loci_required = 0.5,
  alpha = 0.2, max_plot_loci = 500,
  minBin,
  init_m = 0.1,
  num_sweeps = 500,
  burn_in = 100,
  parallel.core = parallel::detectCores() - 1,
  ...) {

  ##test
  # data = "populations.snps.vcf"
  # strata = NULL
  # indivs = "whitelist.samples.tsv"
  # sample.snps = 5000
  # prop_indv_required = 0.5
  # prop_loci_required = 0.5
  # alpha = 0.2
  # max_plot_loci = 5000
  # minBin = 5000
  # init_m = 0.1
  # num_sweeps = 500
  # burn_in = 100
  # parallel.core = parallel::detectCores() - 1
  # random.seed = NULL

  message("whoa: Where's my Heterozygotes? Observations on genotyping Accuracy\n")
  timing.import <- proc.time() # for timing

  if (missing(data)) stop("vcf file missing")
  if (!is.null(strata)) indivs <- NULL # turn off automatically

  # dotslist -------------------------------------------------------------------
  # Note to myself or Eric: If there's a better way to transfer ... from
  # one function to another, I'm all in to lean...


  dotslist <- list(...)
  want <- c("sample.snps", "random.seed")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  whoa.dots <- dotslist[names(dotslist) %in% want]
  sample.snps <- whoa.dots[["sample.snps"]]
  random.seed <- whoa.dots[["random.seed"]]

  # useful outside this function
  if (is.null(random.seed)) {
    random.seed <- sample(x = 1:1000000, size = 1)
    set.seed(random.seed)
  } else {
    set.seed(random.seed)
  }

  v <- read_whoa(
    data = data,
    indivs = indivs,
    filename = NULL,
    parallel.core = parallel.core,
    sample.snps = sample.snps,
    random.seed = random.seed
  )


  if (!is.null(strata)) {
    strata <- readr::read_tsv(file = strata, col_types = c("cc"))
    message("Number of strata to analyze: ", dplyr::n_distinct(strata$STRATA))

    het.miscall <- strata %>%
      split(x = ., f = .$STRATA) %>%
      purrr::map(.x = ., .f = whoa_pop,
                 v = v,
                 indivs = NULL,
                 prop_indv_required = prop_indv_required,
                 prop_loci_required = prop_loci_required,
                 alpha = alpha,
                 max_plot_loci = max_plot_loci,
                 minBin = minBin,
                 init_m = init_m,
                 num_sweeps = num_sweeps,
                 burn_in = burn_in,
                 parallel.core = parallel.core)

    overall.rate <- purrr::map(.x = het.miscall, .f = 3) %>%
      purrr::map(.x = ., .f = 1) %>%
      purrr::map(.x = ., .f = "mean") %>%
      purrr::flatten_dbl(.) %>%
      mean(x = ., na.rm = TRUE) %>%
      round(x = ., digits = 2)

    if (overall.rate > 0.1) {
      message("\n\nOverall heterozygote miscall rate (across read depth and strata): ", overall.rate, " whoa!!!")

    } else {
      message("\n\nOverall heterozygote miscall rate (across read depth and strata): ", overall.rate)
    }

  } else {
    het.miscall <- whoa_pop(
      v = v,
      strata = NULL,
      indivs = NULL,
      prop_indv_required = prop_indv_required,
      prop_loci_required = prop_loci_required,
      alpha = alpha,
      max_plot_loci = max_plot_loci,
      minBin = minBin,
      init_m = init_m,
      num_sweeps = num_sweeps,
      burn_in = burn_in,
      parallel.core = parallel.core)

    overall.rate <- round(x = mean(het.miscall$binned$m_posteriors$mean, na.rm = TRUE), digits = 2)
    if (overall.rate > 0.1) {
      message("\n\nOverall heterozygote miscall rate (across read depth): ", overall.rate, " whoa!!!")
    } else {
      message("\n\nOverall heterozygote miscall rate (across read depth): ", overall.rate)
    }
  }

  message("\nExecution time: ", round(timing.import[[3]]), " sec\n")

  return(res = list(gds = v,
                    het.miscall.rate.analysis = het.miscall,
                    random.seed = random.seed,
                    overall.rate = overall.rate))
  timing.import <- round(proc.time() - timing.import)
}#End run_whoa


# Internal function whoa for 1 pop----------------------------------------------
whoa_pop <- function(
  strata = NULL,
  v,
  indivs = NULL,
  prop_indv_required = 0.5,
  prop_loci_required = 0.5,
  alpha = 0.2,
  max_plot_loci = 500,
  minBin,
  init_m = 0.1,
  num_sweeps = 500,
  burn_in = 100,
  parallel.core = parallel::detectCores() - 1
) {

  # Import VCF -----------------------------------------------------------------
  if (!is.null(strata)) {
    message("\nExecuting whoa for strata: ", unique(strata$STRATA), "...\n")
    indivs <- unique(strata$INDIVIDUALS)
    n.ind <- length(indivs)
    message("Number of samples: ", n.ind)
    vcf.id <- sort(SeqArray::seqGetData(v, "sample.id"))

    # keeping only samples found in VCF...
    blacklisted.id <- setdiff(indivs, vcf.id)
    if (length(blacklisted.id) > 0){
      message("Individuals not in the VCF and excluded from analysis:\n",
              paste(blacklisted.id, collapse = ", "))
    }
    indivs <- intersect(indivs, vcf.id)

    if (length(indivs) == 0) {
      skip.analysis <- TRUE
      message("\nNo sample left to analyse")
      message("Skipping analysis of strata: ", unique(strata$STRATA), "\n")
    } else {
      skip.analysis <- FALSE
    }
    if (!skip.analysis) {
      SeqArray::seqSetFilter(object = v,
                             sample.id = indivs,
                             verbose = FALSE)
    }
  } else {
    skip.analysis <- FALSE
  }
  if (!skip.analysis) {
    # Computing expected and observed genotype frequencies------------------------
    message("Computing expected and observed genotype frequencies")
    gfreqs <- exp_and_obs_geno_freqs(
      v,
      prop_indv_required = prop_indv_required,
      prop_loci_required = prop_loci_required)

    # Plot observed and expected genotype freqs-----------------------------------
    message("Generating expected and observed genotype frequencies plot")
    plot.freqs <- geno_freqs_scatter(gfreqs, alpha = alpha, max_plot_loci = max_plot_loci)
    print(plot.freqs)

    # get posterior estimates for m from different read depth categories ---------
    message("Getting posterior estimates from different read depth categories")

    binned <- infer_m(
      v,
      minBin = minBin,
      indivs = NULL,
      init_m = init_m,
      num_sweeps = num_sweeps,
      burn_in = burn_in)


    # Plot the posterior estimates for heterozyote miscall rates------------------
    message("Generating the plot for the posterior estimates for heterozyote miscall rates")

    plot.post <- posteriors_plot(binned$m_posteriors)
    print(plot.post)

    # reset gds filter -----------------------------------------------------------
    SeqArray::seqResetFilter(v, sample = TRUE, variant = FALSE, verbose = FALSE)
    # we wnat to reset for samples, but not for variant in case a sub-sample was taken...

    # results --------------------------------------------------------------------
    res <- list(gfreqs = gfreqs, plot.freqs = plot.freqs,
                binned = binned, plot.post = plot.post)
  } else {
    res <- list(gfreqs = NULL, plot.freqs = NULL, binned = NULL,
                plot.post = NULL)
  }
  return(res)
}#whoa_pop
