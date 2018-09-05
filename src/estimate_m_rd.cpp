#include <Rcpp.h>
using namespace Rcpp;
#include "whoa.h"

//' Estimate heterozygote miscall rate for different read depth categories (no nulls)
//'
//' To see how this Rcpp function is applied, see the code in
//' \code{\link{infer_m}}.
//' @param Y the 012,-1 matrix that is N x L giving the observed genotypes of the N individuals
//' at L SNPs.
//' @param R integer matrix that is N x L giving the read depth categories.  These must be indexed from
//' 1 up to num_cats.  Missing data should be -1.
//' @param init_m starting value for m.  Typically you might want to use the m estimated
//' from init_m
//' @param num_cats the number of read depth categories.
//' @param p_prior two-vector that holds the beta parameters for a prior on the
//' allele frequency for each locus.  Typical value is c(0.5, 0.5).
//' @param m_prior two-vector that holds the beta parameters for a prior on the
//' heterozygote miscall rate for each locus.  Typical value is c(0.5, 0.5).
//' @param num_reps the number of MCMC sweeps to do.
//' @keywords internal
// [[Rcpp::export]]
List estimate_m_rd(IntegerMatrix Y,
                   IntegerMatrix R,
                   double init_m,
                   int num_cats,
                   NumericVector p_prior,
                   NumericVector m_prior,
                   int num_reps) {

  int i,l, rep;
  int N = Y.nrow();
  int L = Y.ncol();

  // This is to hold the latent genotypes
  IntegerMatrix X(N, L);
  IntegerVector x0(L);   // to count up the number of zero alleles
  IntegerVector x1(L);   // to count up the number of one alleles
  NumericVector p(L); // for reference allele freq
  NumericMatrix Mtrace(num_cats, num_reps + 1);

  // this is to hold the het miscall rates
  NumericVector M(num_cats, init_m);

  // first compute the frequency of the reference allele
  for(l=0;l<L;l++) {
    x0(l) = 0; // this is like having a beta(1/2, 1/2) prior
    x1(l) = 0;
  }


  for(l=0;l<L;l++) {
    for(i=0;i<N;i++) {
      x0(l) += 2 * (Y(i,l) == 0);
      x0(l) += 1 * (Y(i,l) == 1);
      x1(l) += 1 * (Y(i,l) == 1);
      x1(l) += 2 * (Y(i,l) == 2);
    }
    p(l) = (double)(x0(l) + 0.5) / (double)(x0(l) + x1(l) + 1.0);
  }

  // then do the reps, and store the M vector each time
  Mtrace(_, 0) = M;
  for(rep=0; rep<num_reps; rep++) {
    gibbsX(X, Y, R, p, M);
    gibbsP(p, X, p_prior);
    gibbsM(M, num_cats, X, Y, R, m_prior);
    Mtrace(_, rep + 1) = M;
  }


  return List::create(_["M"] = M,
                      _["N"] = N,
                      _["L"] = L,
                      _["p"] = p,
                      _["x0"] = x0,
                      _["x1"] = x1,
                      _["Y"] = Y,
                      _["X"] = X,
                      _["Mtrace"] = Mtrace);
}

