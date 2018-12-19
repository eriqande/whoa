# read vcf for whoa

#' @name read_whoa
#' @title read vcf file for whoa
#' @description read and prepare the vcf file for whoa \href{https://github.com/eriqande/whoa}{whoa}
#' package.


#' @param data (character string) The VCF file. SNPs are biallelic or haplotypes.

#' @param indivs (optional, character) Use this argument to select
#' a subsample of id in the vcf file, e.g. a particular group or population of samples.
#' 2 options: use a character string with desired vcf sample id or a
#' tab-separated file with 1 column, no header, containing the individuals to
#' include in the analysis.
#' By default, this is NULL, in which case everyone from the file is included.
#' Default: \code{indivs = NULL}.

#' @param filename (optional) The file name of the Genomic Data Structure (GDS) file.
#' whoa will append \code{.gds} to the filename.
#' If filename chosen is already present in the
#' working directory, the default \code{whoa_datetime.gds} is chosen.
#' Default: \code{filename = NULL}.

#' @param parallel.core (optional) Default: \code{parallel::detectCores() - 1}.
#' The number of processor used during execution.


#' @param ... (optional) Advance mode that allows to pass further arguments
#' for fine-tuning the function (see details).

#' @return A gds object.

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
#' @rdname read_whoa

#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.
#'
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.


#' @examples
#' \dontrun{
#' require(httr)
#' require(SeqArray)
#' library(whoa)
#'
#' # Get the lobster vcf file:
#' writeBin(httr::content(httr::GET("https://datadryad.org/bitstream/handle/
#' 10255/dryad.108679/10156-586.recode.vcf?sequence=1"), "raw"),"lobster_data_2015.vcf")
#'
#' # Get the strata file containing individuals and groupings
#' writeBin(httr::content(httr::GET("https://datadryad.org/bitstream/handle/
#' 10255/dryad.108679/README.txt?sequence=2"), "raw"),"lobster_strata_2015.tsv")
#'
#' # with built-in defaults:
#' lobster <- whoa::read_whoa(data = "lobster_data_2015.vcf")
#'
#' # Using more arguments: keeping only samples in BUZ sampling sites and using 2000 SNPs:
#' indivs_buz <- readr::read_tsv(file = "lobster_strata_2015.tsv",
#'                           col_names = c("INDIVIDUALS", "STRATA"),
#'                           col_types = "cc") %>%
#'   dplyr::filter(STRATA == "7") %>%
#'   dplyr::select(INDIVIDUALS) %>%
#'   purrr::flatten_chr(.)
#'
#' lobster_buz_2000 <- whoa::read_whoa(
#'     data = "lobster_data_2015.vcf",
#'     indivs = indivs_buz,
#'     sample.snps = 2000,
#'     random.seed = 62994,
#'     filename = "lobster_buz_2000")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


read_whoa <- function(
  data,
  indivs = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  ##Test
  # data = "lobster_data_2015.vcf"
  # indivs = NULL
  # filename <- NULL
  # parallel.core <- parallel::detectCores() - 1
  # sample.snps = 2000
  # random.seed = 62994

  timing.import <- proc.time()
  # Check that SeqArray is installed
  if (!"SeqArray" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install SeqArray for this option:\n
         devtools::install_github("zhengxwen/SeqArray")
         or the bioconductor version:
         source("https://bioconductor.org/biocLite.R")
         biocLite("SeqArray")')
  }

  if (!"gdsfmt" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install gdsfmt for this option:\n
         source("https://bioconductor.org/biocLite.R")
         biocLite("gdsfmt")')
  }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("vcf file missing")
  message("\nReading VCF...")

  # dotslist -------------------------------------------------------------------
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

  if (is.null(random.seed)) {
    random.seed <- sample(x = 1:1000000, size = 1)
    set.seed(random.seed)
  } else {
    set.seed(random.seed)
  }

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    filename <- stringi::stri_join("whoa_", file.date, ".gds")
  } else {
    filename.problem <- file.exists(stringi::stri_join(filename, ".gds"))
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date, ".gds")
    } else {
      filename <- stringi::stri_join(filename, ".gds")
    }
  }

  filename.short <- filename
  filename <- file.path(getwd(), filename)

  # Read vcf -------------------------------------------------------------------
  # Check for bad header generated by software
  check.header <- check_header(data) %$% check.header

  res <- SeqArray::seqVCF2GDS(
    vcf.fn = data,
    out.fn = filename,
    parallel = parallel.core,
    storage.option = "ZIP_RA",
    verbose = FALSE,
    header = check.header,
    info.import = character(0),
    fmt.import = NULL
  ) %>%
    SeqArray::seqOpen(gds.fn = ., readonly = FALSE)

  check.header <- NULL

  gdsfmt::add.gdsn(
    node = res,
    name = "random.seed",
    val = random.seed,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)


  # whitelist samples ----------------------------------------------------------
  if (!is.null(indivs)) {
    if (length(indivs) == 1) {
      indivs <- readr::read_tsv(
        file = indivs, col_names = "INDIVIDUALS", col_types = "c") %>%
        dplyr::arrange(INDIVIDUALS) %>%
        purrr::flatten_chr(.)
    } else {
      indivs <- sort(indivs)
    }

    vcf.id <- sort(SeqArray::seqGetData(res, "sample.id"))
    if (!isTRUE(unique(indivs %in% vcf.id))) {
      stop("samples in indivs don't match the VCF file id's")
    }
    n.ind <- length(indivs)
    SeqArray::seqSetFilter(object = res,
                           sample.id = indivs,
                           verbose = FALSE)
  } else {
    n.ind <- length(SeqArray::seqGetData(res, "sample.id"))
  }

  # sample snps ----------------------------------------------------------------
  if (!is.null(sample.snps)) {
    variant.select <- sample(
      x = SeqArray::seqGetData(res, "variant.id"),
      size = sample.snps)
    SeqArray::seqSetFilter(object = res,
                           variant.id = variant.select,
                           verbose = FALSE)
    n.markers <- length(variant.select)
  } else {
    n.markers <- length(SeqArray::seqGetData(res, "variant.id"))
  }


  # Summary --------------------------------------------------------------------
  message("\nNumber of SNPs: ", n.markers)
  message("Number of samples: ", n.ind)


  message("\nGDS file generated: ", filename.short)
  message("random seed: ", random.seed)
  timing.import <- round(proc.time() - timing.import)
  message("\nWorking time: ", round(timing.import[[3]]), " sec\n")
  return(res)
  } # End read_whoa

