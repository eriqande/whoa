# read vcf for whoa

#' @name read_whoa
#' @title read vcf file into whoa input file
#' @description read vcf file into whoa \href{https://github.com/eriqande/whoa}{whoa}
#' package.


#' @param data (character string) The VCF file. SNPs are biallelic or haplotypes.

#' @param indivs (optional, character) Use this argument to select
#' a subsample of id in the vcf file, e.g. a particular group or population of samples.
#' The file is a tab-separated file with 1 column, no header, containing the individuals to
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

#' \item \code{path.folder}: to write ouput in a specific path
#' (used internally in radiator). Default: \code{path.folder = getwd()}.
#' If the supplied directory doesn't exist, it's created.
#' }

#' @export
#' @rdname read_whoa

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom purrr keep
#' @importFrom readr read_tsv
#' @importFrom stringi stri_join stri_detect_fixed
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
#' # with built-in defaults:
#'  data <- whoa::read_whoa(data = "populations.snps.vcf")
#'
#' # Using more arguments:
#' data <- whoa::read_whoa(
#' data = "populations.snps.vcf",
#' indivs = "select_ind.tsv",
#' parallel.core = 5,
#' filename = "salamander",
#' path.folder = "salamander/het_analysis")
#' }



#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


read_whoa <- function(
  data,
  indivs = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  verbose <- TRUE
  ##Test
  # data = "populations.snps.vcf"
  # filename <- NULL
  # parallel.core <- parallel::detectCores() - 1
  # keep.gds <- TRUE
  # path.folder = NULL

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
  want <- c("keep.gds", "path.folder", "sample.snps", "random.seed")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  radiator.dots <- dotslist[names(dotslist) %in% want]
  keep.gds <- radiator.dots[["keep.gds"]]
  path.folder <- radiator.dots[["path.folder"]]
  sample.snps <- radiator.dots[["sample.snps"]]
  random.seed <- radiator.dots[["random.seed"]]

  # useful outside this function
  if (is.null(keep.gds)) keep.gds <- TRUE
  if (is.null(path.folder)) {
    path.folder <- getwd()
  } else {
    if (!dir.exists(path.folder)) dir.create(path.folder)
  }

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
  filename <- file.path(path.folder, filename)

  # Read vcf -------------------------------------------------------------------
  # Check for bad header generated by software
  detect.source <- check_header(data)
  check.header <- detect.source$check.header

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

  check.header <- detect.source <- NULL

  gdsfmt::add.gdsn(
    node = res,
    name = "random.seed",
    val = random.seed,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)


  # whitelist samples ----------------------------------------------------------
  if (!is.null(indivs)) {
    indivs <- readr::read_tsv(
      file = indivs, col_names = "INDIVIDUALS", col_types = "c") %>%
      dplyr::arrange(INDIVIDUALS)
    vcf.id <- sort(SeqArray::seqGetData(res, "sample.id"))
    if (!isTRUE(unique(indivs$INDIVIDUALS %in% vcf.id))) {
      stop("samples in the whitelist don't match the VCF file id's")
    }
    n.ind <- length(indivs$INDIVIDUALS)
    SeqArray::seqSetFilter(object = res,
                           sample.id = indivs$INDIVIDUALS,
                           verbose = FALSE)
  } else {
    n.ind <- length(SeqArray::seqGetData(res, "sample.id"))
  }

  # sample snps ----------------------------------------------------------------
  if (!is.null(sample.snps)) {
    message("random seed: ", random.seed)
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

  if (verbose && keep.gds) {
    message("\nGDS file generated: \n", filename.short)
  }
  timing.import <- round(proc.time() - timing.import)
  if (verbose) message("\nWorking time: ", round(timing.import[[3]]), " sec\n")
  return(res)
  } # End read_whoa

# Internal nested Function -----------------------------------------------------
# check_header
#' @title Check the vcf header and detect vcf source
#' @description Check the vcf header and detect vcf source
#' @rdname check_header
#' @keywords internal
#' @export
check_header <- function(vcf) {
  check.header <- SeqArray::seqVCF_Header(vcf)

  if (check.header$format$Number[check.header$format$ID == "AD"] == 1) {
    check.header$format$Number[check.header$format$ID == "AD"] <- "."
  }

  problematic.id <- c("AD", "AO", "QA", "GL")
  problematic.id <- purrr::keep(.x = problematic.id, .p = problematic.id %in% check.header$format$ID)
  for (p in problematic.id) {
    check.header$format[check.header$format$ID == p, "Number"] <- "."
  }
  # check.header$format

  check.source <- check.header$header$value[check.header$header$id == "source"]
  is.stacks <- stringi::stri_detect_fixed(str = check.source, pattern = "Stacks")
  if (is.stacks) {
    stacks.2 <- keep.stacks.gl <- stringi::stri_detect_fixed(
      str = check.source,
      pattern = "Stacks v2")
    if (!keep.stacks.gl) {
      check.header$format <- dplyr::filter(check.header$format, ID != "GL")
    }
      } else {
    stacks.2 <- FALSE
  }
  return(res = list(source = stacks.2, check.header = check.header))
}#End check_header
