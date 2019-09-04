

#' Restriction-associated digest (RAD) sequence data from 36 lobsters at 2000 single nucleotide polymorphisms (SNPs)
#'
#' A data set from the study by Benestan et al. (2016).
#'
#' @format A `vcfR' object in which the original data set has
#' been reduced to just the 36 lobsters from the BUZ population
#' and a randomly sampled 2000 SNPs from the 10,156 originally
#' available.
#'
#' This is the sort of object obtained after calling \code{vcfR::read.vcfR()}
#' on a variant call format (VCF) file.
#'
#' @source \url{https://datadryad.org//resource/doi:10.5061/dryad.q771r.3}
"lobster_buz_2000"


#' An 012 matrix from  (RAD) sequence data from 36 lobsters at 2000 single nucleotide polymorphisms (SNPs)
#'
#' A data set from the study by Benestan et al. (2016).
#'
#' @format This is an integer matrix with positions in rows and samples in columns.
#' It is an 012 matrix that corresponds to lobster_buz_2000.  It is useful as
#' an example of the necessary format for the d012 argument to \code{\link{exp_and_obs_geno_freqs}}.
#'
#'
#' @source \url{https://datadryad.org//resource/doi:10.5061/dryad.q771r.3}
"lobster_buz_2000_as_012_matrix"
