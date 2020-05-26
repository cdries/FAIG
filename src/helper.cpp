#include "RcppArmadillo.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @export get_valmat
// [[Rcpp::export]]
arma::mat get_valmat(arma::mat vals, arma::ivec alloc, int n_items, int n_persons) {
  // computes the matrix of dimension n_person x n_persons where each row contains
  // the valuation of that person for the different sets of items the other persons
  // receive
  // 
  // arguments:
  // vals     : matrix (n_persons x n_items) with each row the valuation of that person for the items
  // alloc    : index of the person to which each item belongs (1, 2, ..., n_persons)
  // n_items  : number of items (length of alloc)
  // n_persons : number of persons
  //
  // output:
  // valmat   : valuation matrix of the different sets (columns) to each person (row)
  //
  // author: Dries Cornilly

  arma::mat valmat = arma::zeros(n_persons, n_persons);
  for (int ii = 0; ii < n_items; ii++) {
    valmat.col(alloc(ii) - 1) += vals.col(ii);
  }
  
  return valmat;
}


//' @export
// [[Rcpp::export]]
double get_maxenvy(arma::mat valmat, int n) {
  // gets the maximum envy from a matrix with valuations
  // 
  // arguments:
  // valmat   : valuation matrix of the different sets (columns) to each person (row)
  // n        : number of rows / columns of valmat
  //
  // output:
  // maxenvy  : maximum envy of the allocation
  //
  // author: Dries Cornilly
  
  for (int ii = 0; ii < n; ii++) {
    valmat.row(ii) -= valmat(ii, ii);
  }
  double maxenvy = valmat.max();
  
  return maxenvy;
}


//' @export
//[[Rcpp::export]]
arma::mat get_avgval(arma::mat valmat, int n) {
  // gets the average value of the valuations of each person, arranged to have the same
  // dimensions as valmat - valmat can also be vals matrix of dimension n_persons x n_items,
  // in this case, n is still the number of persons
  // 
  // arguments:
  // valmat   : valuation matrix of the different sets (columns) to each person (row)
  // n        : number of rows / columns of valmat - number of persons
  //
  // output:
  // avgval   : matrix with average valuations for each agent (row i: agent i, repeat n times)
  //
  // author: Dries Cornilly
  
  double nn = n;
  arma::mat avgval = arma::repmat(arma::sum(valmat, 1) / nn, 1, n);
  return avgval;
}


//' @export
//[[Rcpp::export]]
double get_fnV(arma::mat valmat, int n, arma::mat avgval) {
  // gets the social inequality measure from a matrix with valuations
  // 
  // arguments:
  // valmat   : valuation matrix of the different sets (columns) to each person (row)
  // n        : number of rows / columns of valmat
  // avgval   : matrix with average valuations for each person (all columns are the same as the first)
  //
  // output:
  // v        : social inequality measure
  //
  // author: Dries Cornilly
  
  double nn = 1.0 * n * n;
  double v = arma::sum(arma::sum(arma::square(valmat - avgval), 1)) / nn;
  return v;
}
