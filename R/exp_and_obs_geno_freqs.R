

#' Computed expected and observed genotype frequencies from a `vcfR' object
#'
#' Under the assumption of Hardy-Weinberg equilibrium, this function uses the
#' estimated allele frequencies from the data set in v to compute expected genotype
#' frequencies, and then reports these along with the observed genotype frequencies.
#' Loci come out named as CHROM--POS.
#'
#' @param v a `vcfR' object.  Exactly one of \code{v} or \code{d012} is required. But you
#' can't use both!
#' @param d012 an integer matrix (or a numeric matrix, which will be coerced to
#' be of integer type) with individuals in columns, and markers in rows.
#' 0 denotes a genotype homozygous for the reference allele, 1 is a heterozygote, 2 is a
#' homozygote for the alternate allele, and -1 denotes missing data.  This matrix is not
#' required to have column (sample) names.  They won't be used if they are present.  But,
#' the matrix must have rownames, which should be in the format of CHROM--POS (i.e. the "chromosome"
#' name (or the "contig" name) followed by a "--" followed by the position of the marker in the "chromosome").
#' Exactly one of \code{v} or \code{d012} is required. But you
#' can't use both!
#' @param prop_indv_required loci will be dropped if a proportion of
#' individuals less than prop_indv_required have non-missing data at that locus.
#' Default is 0.5
#' @param prop_loci_required individual will be dropped if their proportion of
#' non-missing loci is less than prop_loci_required. Default is 0.5.
#' @export
#' @return Returns a tibble with the following columns: \code{snp} = the locus name
#' as CHROM--POS; \code{p} = The frequency of the alternate (ALT) allele; \code{ntot} = the total
#' number of individuals with no missing data at the locus; \code{geno} = column
#' telling which genotype (0, 1, or 2) is referred to; \code{p_exp} = expected
#' frequency of the genotype; \code{p_obs} = observed frequency of genotype;
#' \code{n_exp} = expected number of such genotypes; \code{n_obs} = observed
#' number of such genotypes; \code{z_score} = simple statistic giving how far
#' the observed genotype frequency is from that expected under Hardy-Weinberg
#' equilibrium.
#'
#' @examples
#' eao <- exp_and_obs_geno_freqs(v = lobster_buz_2000)
#'
#' # if you wanted to run that on an 012 matrix,
#' # it would be like this:
#' eao012 <- exp_and_obs_geno_freqs(d012 = lobster_buz_2000_as_012_matrix)
exp_and_obs_geno_freqs <- function(v = NULL, d012 = NULL, prop_indv_required = 0.5, prop_loci_required = 0.5) {

  if(!xor(is.null(v), is.null(d012))) {
    stop("Exactly one of v or d012 must be present.")
  }

  if(!is.null(v)) {
    # first get an 012 matrix with loci in the columns
    # and samples in the rows
    tmp <- vcfR::extract.gt(v, element = "GT")
    dimnames(tmp) <- NULL
    g012 <- t(make_it_012(tmp))
    colnames(g012) <- paste(v@fix[,1], v@fix[,2], sep = "--")  # name the Chrom--Pos
  }
  if(!is.null(d012)) {
    # check that it is numeric or integer
    if(!(storage.mode(d012) == "integer" | storage.mode(d012) == "numeric")) {
      stop("d012 must be a matrix of type integer or numeric, not ", storage.mode(d012))
    }
    # make sure to set storage mode to integer
    storage.mode(d012) <- "integer"

    # check to make sure there are colnames
    if(is.null(rownames(d012))) {
      stop("You must have rownames on d012")
    }
    # remove the colnames if they are there
    colnames(d012) <- NULL

    # then return g012 as the transpose of d012
    g012 <- t(d012)
  }

  # now, drop the loci with too much missing data, and after that
  # drop the indivs with too much missing data
  pml <- colMeans(g012 > -1)
  g012 <- g012[, pml > prop_indv_required]
  pmi <- rowMeans(g012 > -1)
  g012 <- g012[pmi > prop_loci_required,]

  gfc <- geno_freq_calc_single(g012)

  gfc
}
