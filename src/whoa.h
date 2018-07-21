#include <Rcpp.h>
using namespace Rcpp;


// some function prototypes
void gibbsM(NumericVector M, int num_cats, IntegerMatrix X, IntegerMatrix Y, IntegerMatrix R, NumericVector pri);
void gibbsP(NumericVector p, IntegerMatrix X, NumericVector pri);
void gibbsX(IntegerMatrix X, IntegerMatrix Y, IntegerMatrix R, NumericVector p, NumericVector M);
