#include "RcppArmadillo.h"
#include "helper.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List mincov(arma::mat vals, arma::ivec alloc, arma::mat beta, int maxiter, int maxnoimprove, double eps) {
  // mincov algorithm - randomly choose a column (item) and give it to the person (row) that should 
  // receive it according to the theorem in the paper. This is done a maximum of maxiter steps, 
  // until a social inequality of eps, or until there is no improvement for maxnoimprove steps.
  // 
  // arguments:
  // vals     : matrix (n_persons x n_items) with each row the valuation of that person for the items
  // alloc    : index of the person to which each item belongs
  // beta     : beta of each person and item with respect to the first person
  // maxiter  : maximum number of iterations
  // maxnoimprove : terminate if no improvement for maxnoimprove consecutive steps
  // eps      : terminate if maxenvy < eps
  //
  // output:
  // alloc    : optimal allocation
  // minsoc   : optimal value of social inequality - corresponds to alloc
  // valmat   : valuation matrix of the different sets (columns) to each person (row)
  // socvec   : vector with social inequality values at the different iterations
  // status   : 0 (V < eps); 1 (maxiter reached)
  // iter     : number of iterations the algorithm completed before stopping
  //
  // author: Dries Cornilly
  
  // initialize
  int n_items = vals.n_cols;                      // number of items to distribute
  int n_persons = vals.n_rows;                    // number of persons to distribute among
  arma::mat valmat = get_valmat(vals, alloc, n_items, n_persons); // get value of each set of items for each person
  arma::mat avgval = get_avgval(vals, n_persons); // initialize average valuations
  arma::vec socvec = arma::zeros(maxiter + 1);    // initialize social inequality through iterations
  socvec(0) = get_fnV(valmat, n_persons, avgval);

  // iterate
  int iter = 0;
  bool converged = false;
  int status = 1;
  int noimprove = 0;
  while (iter < maxiter && !converged) {

    // sample item and get its owner
    int item = arma::randi(1, arma::distr_param(0, n_items - 1))(0);
    int oldperson = alloc(item) - 1;

    // determine whom to give it to
    arma::mat valmattemp = valmat;
    valmattemp.col(oldperson) -= vals.col(item);
    for (int jj = 0; jj < n_persons; jj++) valmattemp.row(jj) *= beta(jj, item);
    arma::rowvec L = arma::sum(valmattemp, 0);
    int newperson = arma::index_min(L);

    // give item
    valmat.col(oldperson) -= vals.col(item);
    valmat.col(newperson) += vals.col(item);
    alloc(item) = newperson + 1;
    socvec(1 + iter) = get_fnV(valmat, n_persons, avgval);
    if (socvec(1 + iter) < socvec(iter)) {
      noimprove = 0;
    } else {
      noimprove++;
    }

    // check convergence
    if (socvec(1 + iter) < eps) {
      converged = true;
      status = 0;
    } else if (noimprove >= maxnoimprove) {
      converged = true;
      status = 2;
    }

    iter++;
  }

  List out;
  out["alloc"] = alloc;
  out["minsoc"] = socvec(iter);
  out["valmat"] = valmat;
  out["socvec"] = socvec;
  out["status"] = status;
  out["iter"] = iter;
  
  return out;
}


// [[Rcpp::export]]
List mincovtarget(arma::mat vals, arma::ivec alloc, arma::mat beta, arma::vec target, 
                  int maxiter, int maxnoimprove, double eps) {
  // mincov algorithm with target value - randomly choose a column (item) and give it to the person (row) 
  // that should receive it according to the theorem in the paper. This is done a maximum of maxiter steps, 
  // until a social inequality of eps, or until there is no improvement for maxnoimprove steps.
  // 
  // arguments:
  // vals     : matrix (n_persons x n_items) with each row the valuation of that person for the items
  // alloc    : index of the person to which each item belongs
  // beta     : beta of each person and item with respect to the first person
  // target   : target value for each of the persons
  // maxiter  : maximum number of iterations
  // maxnoimprove : terminate if no improvement for maxnoimprove consecutive steps
  // eps      : terminate if maxenvy < eps
  //
  // output:
  // alloc    : optimal allocation
  // minsoc   : optimal value of social inequality - corresponds to alloc
  // valmat   : valuation matrix of the different sets (columns) to each person (row)
  // socvec   : vector with social inequality values at the different iterations
  // status   : 0 (V < eps); 1 (maxiter reached)
  // iter     : number of iterations the algorithm completed before stopping
  //
  // author: Dries Cornilly
  
  // initialize
  int n_items = vals.n_cols;                      // number of items to distribute
  int n_persons = vals.n_rows;                    // number of persons to distribute among
  arma::mat valmatT = get_valmat(vals, alloc, n_items, n_persons); // get value of each set of items for each person
  valmatT.diag() -= target;
  arma::mat avgval = get_avgval(valmatT, n_persons); // initialize average valuations
  arma::vec socvec = arma::zeros(maxiter + 1);    // initialize social inequality through iterations
  socvec(0) = get_fnV(valmatT, n_persons, avgval);

  // iterate
  int iter = 0;
  bool converged = false;
  int status = 1;
  int noimprove = 0;
  while (iter < maxiter && !converged) {
    
    // sample item and get its owner
    int item = arma::randi(1, arma::distr_param(0, n_items - 1))(0);
    int oldperson = alloc(item) - 1;
    
    // determine whom to give it to
    arma::mat valmattemp = valmatT;
    valmattemp.col(oldperson) -= vals.col(item);
    for (int jj = 0; jj < n_persons; jj++) valmattemp.row(jj) *= beta(jj, item);
    arma::rowvec L = arma::sum(valmattemp, 0);
    int newperson = arma::index_min(L);
    
    // give item
    valmatT.col(oldperson) -= vals.col(item);
    valmatT.col(newperson) += vals.col(item);
    alloc(item) = newperson + 1;
    socvec(1 + iter) = get_fnV(valmatT, n_persons, avgval);
    if (socvec(1 + iter) < socvec(iter)) {
      noimprove = 0;
    } else {
      noimprove++;
    }
    
    // check convergence
    if (socvec(1 + iter) < eps) {
      converged = true;
      status = 0;
    } else if (noimprove >= maxnoimprove) {
      converged = true;
      status = 2;
    }

    iter++;
  }
  
  // compute end-statistics without the target columns
  arma::mat valmat0 = get_valmat(vals, alloc, n_items, n_persons);
  arma::mat avgval0 = get_avgval(valmat0, n_persons);
  double minsoc0 = get_fnV(valmat0, n_persons, avgval0);
  
  List out;
  out["alloc"] = alloc;
  out["minsocT"] = socvec(iter);
  out["minsoc"] = minsoc0;
  out["valmatT"] = valmatT;
  out["valmat"] = valmat0;
  out["socvec"] = socvec;
  out["status"] = status;
  out["iter"] = iter;
  
  return out;
}
