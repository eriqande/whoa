% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_depths.R
\name{bin_depths}
\alias{bin_depths}
\title{bin read depths of SNPs into categories having at least S observations}
\usage{
bin_depths(D, S)
}
\arguments{
\item{D}{a matrix of read depths.  Rows are individuals, columns are SNPs.  Cells where data are missing
in the genotype matrix must be denoted as NA}

\item{S}{the min number of observations to have in each bin}
}
\value{
This returns a list with two components.  \code{dp_bins} is a matrix of the same
shape as D with the bin categories (as 1, 2, ...) and -1 for this cells
corresponding to missing genotypes.  \code{num_cats} is the number of depth bins.
\code{tidy_bins} is a long format description of the bins.
\code{bin_stats} is a tibble giving summary information about the read depth bins which
is useful for plotting things, etc.
}
\description{
bin read depths of SNPs into categories having at least S observations
}
\examples{

# get a matrix of read depths and make it an integer matrix
depths <- vcfR::extract.gt(lobster_buz_2000, element = "DP")
storage.mode(depths) <- "integer"

# get a character matrix of genotypes, so we can figure out which
# are missing and mask those from depths
genos <- vcfR::extract.gt(lobster_buz_2000, element = "GT")

# make missing in depths if missing in genos
depths[is.na(genos)] <- NA

# bin the read depths into bins with at least 1000 observations in each bin
bins <- bin_depths(depths, 1000)
}
\keyword{internal}
