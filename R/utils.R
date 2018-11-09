# Global Variables -------------------------------------------------------------

# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if(getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "Freq",
      "dp",
      "bin",
      "mean_dp",
      "lo95",
      "hi95",
      "total_n",
      "0",
      "1",
      "2",
      "geno",
      "n_exp",
      "n_obs",
      "ntot",
      "p_exp",
      "p_obs",
      "snp",
      "z_score",
      "INDIVIDUALS",
      "ID"
    )
  )
}

# magrittr ---------------------------------------------------------------------
#' @title Forward-pipe operator
#' @description magrittr forward-pipe operator
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

# Exposition pipe-operator
#' @title Exposition pipe-operator
#' @description magrittr Exposition pipe-operator
#' @name %$%
#' @rdname Exposition_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %$%
#' @usage lhs \%$\% rhs
NULL

# compound assignment pipe operator
#' @title compound assignment pipe operator
#' @description magrittr compound assignment pipe operator
#' @name %<>%
#' @rdname compound_assignment_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
NULL

# check_header internal function in read_whoa.R --------------------------------
#' @title Check the vcf header and detect vcf source
#' @description Check the vcf header and detect vcf source
#' @rdname check_header
#' @keywords internal
#' @export
check_header <- function(vcf) {
  check.header <- SeqArray::seqVCF_Header(vcf)
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


