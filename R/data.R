

#' RAD-seq data from 36 lobsters at 2000 SNPs
#'
#' A data set from the study by Benestan et al. (2016).
#'
#' @format A vcfR object in which the original data set has
#' been reduced to just the 36 lobsters from the BUZ population
#' and a randomly sampled 2000 SNPs from the 10,156 originally
#' available.
#'
#' This is the sort of object obtained after calling \code{vcfR::read.vcfR()}
#' on a VCF file.
#' @references Benestan L, Gosselin T, Perrier C et al. (2015)
#' RAD genotyping reveals fine-scale genetic structuring and provides powerful
#' population assignment in a widely distributed marine species,
#' the American lobster (Homarus americanus). Molecular Ecology, 24, 3299-3315.
#'
#' @source \url{https://datadryad.org//resource/doi:10.5061/dryad.q771r.3}
"lobster_buz_2000"
