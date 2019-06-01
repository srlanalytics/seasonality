// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec getrow(arma::field<rowvec> A){
  rowvec b = A(2);
  uword elm = A.n_elem;
  Rcpp::Rcout << elm << endl;
  // List out;
  // out["b"] = b;
  // out["elm"] = elm;
  return(b);
}
