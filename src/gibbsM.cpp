#include <Rcpp.h>
using namespace Rcpp;

//' simulate a new miscall rate for each read depth category given X and Y
//'
//' This just writes new values into M as if it were an output variable
//' @keywords internal
// [[Rcpp::export]]
void gibbsM(NumericVector M, int num_cats, IntegerMatrix X, IntegerMatrix Y, IntegerMatrix R, NumericVector pri) {
  int i,l,k;
  int N = Y.nrow();
  int L = Y.ncol();
  int r, x, y;

  NumericVector mc(num_cats, pri(0));  // for counting up the number of hets that were miscalled
  NumericVector cc(num_cats, pri(1));  // for counting up the number of hets that were correctly called

  // Note, the above two vectors are initialized to the beta prior for the miscall rate

  // cycle over all loci and all individuals
  for(l=0;l<L;l++) {
    for(i=0;i<N;i++) {
      r = R(i, l) - 1;  // the read depth category
      y = Y(i, l);
      x = X(i, l);
      if(r >= 0 && y >= 0) { // as long as genotype is not missing and read depth category is valid
        if(x == 1 && y != 1) mc(r) += 1.0;
        if(x == 1 && y == 1) cc(r) += 1.0;
      }
    }
  }

  // now simulate new values from a beta
  for(k=0;k<num_cats;k++) {
    M(k) = R::rbeta(mc(k), cc(k));
  }

}
