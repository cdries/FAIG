#ifndef HELPER_H
#define HELPER_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::mat get_valmat(arma::mat vals, arma::ivec alloc, int n_items, int n_persons);

double get_maxenvy(arma::mat valmat, int n);

arma::mat get_avgval(arma::mat valmat, int n);

double get_fnV(arma::mat valmat, int n, arma::mat avgval);

double get_util(arma::mat valmat);


#endif
