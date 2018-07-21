#include <Rcpp.h>
using namespace Rcpp;




//' compute full conditional for each X given Y, p, R, and m, and then sample from it
//'
//' This just writes new values into X as if it were an output variable
//' @keywords internal
// [[Rcpp::export]]
void gibbsX(IntegerMatrix X, IntegerMatrix Y, IntegerMatrix R, NumericVector p, NumericVector M) {
  int i,l,k;
  int N = Y.nrow();
  int L = Y.ncol();
  double probs[3];  // for storing the probs
  double pp, m, normo, rr, cumul;
  int y, r, newX;

  // cycle over each individual and locus
  for(i=0;i<N;i++) {
    for(l=0;l<L;l++) {
      pp = p(l);  // get the reference alle freq
      // intialize the probs to the prior from the allele freq
      probs[0] = pp * pp;
      probs[1] = 2.0 * pp * (1 - pp);
      probs[2] = (1 - pp) * (1 - pp);

      y = Y(i,l); // temp variable to hold genotype
      r = R(i,l) - 1;   // temp variable to hold read category.  Subtract 1 for base-0 indexing...
      if(y >= 0 && r >= 0) {  // if we have a genotype call and a read depth category then use those to update the priors
        m = (double)(M(r));
        probs[0] *= (y == 0);
        probs[1] *= ( (m / 2.0) * (y == 0) + (1.0 - m) * (y == 1) + (m / 2.0) * (y == 2));
        probs[2] *= (y == 2);
      }

      // now normalize those
      normo = 0.0;
      for(k=0;k<3;k++) normo += probs[k];
      for(k=0;k<3;k++) probs[k] /= normo;

      // now sample from them
      rr = R::runif(0,1);;
      cumul = 0.0;
      newX = -2;
      for(k=0;k<3;k++) {
        cumul += probs[k];
        if(cumul > rr) {
          newX = k;
          break;
        }
      }

      // and assign that to X(i,l)
      X(i,l) = newX;

    }
  }

}
