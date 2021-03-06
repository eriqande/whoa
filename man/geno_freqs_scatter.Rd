% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geno_freqs_scatter.R
\name{geno_freqs_scatter}
\alias{geno_freqs_scatter}
\title{return a `ggplot2' plot object of observed and expected genotype freqs}
\usage{
geno_freqs_scatter(gfc, alpha = 0.2, max_plot_loci = 500)
}
\arguments{
\item{gfc}{a tibble like that created by exp_and_obs_geno_freqs()}

\item{alpha}{the transparency (alpha) parameter to apply to the points
in the scatterplot. Default is 0.2.}

\item{max_plot_loci}{By default this plots only 500 loci, sampled
randomly, to keep `ggplot2' taking forever to plot, for example, 100K
points.  If you want to plot all the points, set this to a number
larger than the number of single nucleotide polymorphisms (SNPs) in the data set.}
}
\description{
return a `ggplot2' plot object of observed and expected genotype freqs
}
\examples{
# get the expected and observed geno freqs
gfreqs <- exp_and_obs_geno_freqs(lobster_buz_2000)
g <- geno_freqs_scatter(gfreqs)

# now g is a 'ggplot2' object.
}
