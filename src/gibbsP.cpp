#include <Rcpp.h>
using namespace Rcpp;


//' simulate new reference allele frequencies from their beta full conditional
//'
//' This just writes new values into P as if it were an output variable
//' @keywords internal
// [[Rcpp::export]]
void gibbsP(NumericVector p, IntegerMatrix X, NumericVector pri) {
  int i,l;
  int N = X.nrow();
  int L = X.ncol();

  double x0;   // to count up the number of zero alleles
  double x1;   // to count up the number of one alleles

  for(l=0;l<L;l++) {
    x0 = pri(0);  // initialize to the prior values
    x1 = pri(1);
    for(i=0;i<N;i++) {
      x0 += 2.0 * (X(i,l) == 0) + 1.0 * (X(i,l) == 1);
      x1 += 2.0 * (X(i,l) == 2) + 1.0 * (X(i,l) == 1);
    }
    p(l) = R::rbeta(x0, x1);
  }
}

