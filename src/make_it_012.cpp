#include <Rcpp.h>
using namespace Rcpp;



//' just a quick function for making an 012 matrix from a character matrix
//'
//' The standard way within R of pulling values out of a named
//' vector really bogs down on large data sets.  So I will do this instead.
//' @param M a character matrix of VCF genotypes and no dimnames.  Allowable values are
//' "0/0", and "0|0", which get coverted to integer 0;  "0/1", "0|1", "1/0", and "1|0",
//' which get converted to integer 1; and
//'  "1/1", and "1|1", which get converted to integer 2.  Everything else gets
//'  converted to -1 to denote missing data.
//' @return An integer matrix of values which are 0, 1, 2, or -1.
//' @export
//' @examples
//' # get an 012 matrix from the lobster data
//' tmp <- t(vcfR::extract.gt(lobster_buz_2000, element = "GT"))
//' dimnames(tmp) <- NULL
//' g <- make_it_012(tmp)
// [[Rcpp::export]]
IntegerMatrix make_it_012(CharacterMatrix M) {
  int N = M.ncol();
  int L = M.nrow();
  int i,l;

  IntegerMatrix ret(L, N);

  for(l=0;l<L;l++) {
    for(i=0;i<N;i++) {
      if(strcmp(M(l,i), "0/0") == 0) {
        ret(l,i) = 0L;
      } else if(strcmp(M(l,i), "0|0") == 0) {
        ret(l,i) = 0L;
      } else if(strcmp(M(l,i), "0/1") == 0) {
        ret(l,i) = 1L;
      } else if(strcmp(M(l,i), "1/0") == 0) {
        ret(l,i) = 1L;
      } else if(strcmp(M(l,i), "0|1") == 0) {
        ret(l,i) = 1L;
      } else if(strcmp(M(l,i), "1|0") == 0) {
        ret(l,i) = 1L;
      } else if(strcmp(M(l,i), "1/1") == 0) {
        ret(l,i) = 2L;
      } else if(strcmp(M(l,i), "1|1") == 0) {
        ret(l,i) = 2L;
      } else {
        ret(l,i) = -1L;  // this takes care of NAs
      }
    }
  }

  return(ret);
}









