#include "RcppArmadillo.h"
#include "helper.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::ivec random_alloc(int n_items, int n_persons) {
  // generate random allocation of the item to the persons
  //
  // arguments:
  // n_items  : number of items (length of alloc)
  // n_persons : number of persons
  //
  // output:
  // alloc    : index of the person to which each item belongs (in 1, 2, ..., n_persons)
  //
  // author: Dries Cornilly
  
  arma::ivec alloc = arma::randi(n_items, arma::distr_param(1, n_persons));
  return alloc;
}


// [[Rcpp::export]]
List randselect_envy(arma::mat vals, int maxiter, double eps) {
  // random minmaxenvy algorithm - randomly (uniformly) allocate each item to one of the persons. This
  // is done a maximum of maxiter steps or until maxenvy of eps is reached.
  // 
  // arguments:
  // vals     : matrix (n_persons x n_items) with each row the valuation of that person for the items
  // maxiter  : maximum number of iterations
  // eps      : terminate if maxenvy < eps
  //
  // output:
  // alloc    : optimal allocation 
  // minmaxenvy : optimal value of maxenvy - corresponds to alloc
  // valmat   : valuation matrix of the different sets (columns) to each person (row)
  // envyvec  : vector with maxenvy values at the different iterations
  // status   : 0 (V < eps); 1 (maxiter reached)
  // iter     : number of iterations the algorithm completed before stopping
  //
  // author: Dries Cornilly
  
  // initialize
  int n_items = vals.n_cols;                      // number of items to distribute
  int n_persons = vals.n_rows;                    // number of persons to distribute among
  double minmaxenvy = arma::sum(arma::sum(vals)); // initialize at high value
  arma::ivec alloc(n_items);
  arma::mat valmat(n_persons, n_persons);
  arma::vec envyvec = arma::zeros(maxiter);
  
  // iterate
  int iter = 0;
  bool converged = false;
  int status = 1;
  while (iter < maxiter && !converged) {
    
    // give items to random person
    arma::ivec alloctemp = random_alloc(n_items, n_persons);
    arma::mat valmattemp = get_valmat(vals, alloctemp, n_items, n_persons);
    double envytemp = get_maxenvy(valmattemp, n_persons);
    envyvec(iter) = envytemp;
    
    // update optimal solution
    if (envytemp < minmaxenvy) {
      minmaxenvy = envytemp;
      alloc = alloctemp;
      valmat = valmattemp;
    }
    
    // check convergence
    if (minmaxenvy < eps) {
      converged = true;
      status = 0;
    }
    
    iter++;
  }
  
  List out;
  out["alloc"] = alloc;
  out["minmaxenvy"] = minmaxenvy;
  out["valmat"] = valmat;
  out["envyvec"] = envyvec;
  out["status"] = status;
  out["iter"] = iter;
  
  return out;
}


// [[Rcpp::export]]
List randselect_social(arma::mat vals, int maxiter, double eps) {
  // random social inequality algorithm - randomly (uniformly) allocate each item to one of the
  // persons. This is done a maximum of maxiter steps or until a social inequality of eps is reached.
  //
  // arguments:
  // vals     : matrix (n_persons x n_items) with each row the valuation of that person for the items
  // maxiter  : maximum number of iterations
  // eps      : terminate if soc < eps
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
  double minsoc = arma::sum(arma::sum(vals % vals)); // initialize at high value
  arma::ivec alloc(n_items);
  arma::mat avgval = get_avgval(vals, n_persons);
  arma::mat valmat(n_persons, n_persons);
  arma::vec socvec = arma::zeros(maxiter);

  // iterate
  int iter = 0;
  bool converged = false;
  int status = 1;
  while (iter < maxiter && !converged) {

    // give items to random person
    arma::ivec alloctemp = random_alloc(n_items, n_persons);
    arma::mat valmattemp = get_valmat(vals, alloctemp, n_items, n_persons);
    double soctemp = get_fnV(valmattemp, n_persons, avgval);
    socvec(iter) = soctemp;

    // update optimal solution
    if (soctemp < minsoc) {
      minsoc = soctemp;
      alloc = alloctemp;
      valmat = valmattemp;
    }

    // check convergence
    if (minsoc < eps) {
      converged = true;
      status = 0;
    }

    iter++;
  }

  List out;
  out["alloc"] = alloc;
  out["minsoc"] = minsoc;
  out["valmat"] = valmat;
  out["socvec"] = socvec;
  out["status"] = status;
  out["iter"] = iter;

  return out;
}
