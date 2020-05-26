#include "RcppArmadillo.h"
#include "helper.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
int testfunc(int oldperson, int addperson, int n_persons){
  int newperson = (oldperson + addperson) % n_persons;
  return newperson;
}

// [[Rcpp::export]]
List localtrades_envy(arma::mat vals, arma::ivec alloc, int maxiter, double eps) {
  // envy-swapping algorithm - randomly choose an item and allocate it to a different person, if it
  // decreases the maxenvy objective. This is done a maximum of maxiter steps or until a maxenvy of eps is
  // reached.
  // 
  // arguments:
  // vals     : matrix (n_persons x n_items) with each row the valuation of that person for the items
  // alloc    : index of the person to which each item belongs
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
  arma::mat valmat = get_valmat(vals, alloc, n_items, n_persons); // get value of each set of items for each person
  arma::vec envyvec = arma::zeros(maxiter + 1);   // initialize maxenvy through iterations
  double minmaxenvy = get_maxenvy(valmat, n_persons); // maxenvy at initial stage
  envyvec(0) = minmaxenvy;
  
  // iterate
  int iter = 0;
  bool converged = false;
  int status = 1;
  while (iter < maxiter && !converged) {
    
    // sample items to give to a different owner
    int item = arma::randi(1, arma::distr_param(0, n_items - 1))(0);
    int addperson = arma::randi(1, arma::distr_param(1, n_persons - 1))(0);
    
    // try the reassignment and update if improvements are made
    arma::ivec alloctemp = alloc - 1;
    arma::mat valmattemp = valmat;
    int oldperson = alloctemp(item);
    int newperson = (oldperson + addperson) % n_persons;
    valmattemp.col(oldperson) -= vals.col(item);
    valmattemp.col(newperson) += vals.col(item);
    alloctemp(item) = newperson;
    double envytemp = get_maxenvy(valmattemp, n_persons);
    envyvec(1 + iter) = envytemp;
    
    // update if lower maxenvy
    if (envytemp < minmaxenvy) {
      minmaxenvy = envytemp;
      alloc = alloctemp + 1;
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
List localtrades_social(arma::mat vals, arma::ivec alloc, int maxiter, double eps) {
  // social inequality-swapping algorithm - randomly choose an item and allocate it to a different 
  // person, if it decreases the social inequality objective. This is done a maximum of maxiter steps 
  // or until a maxenvy of eps is reached.
  // 
  // arguments:
  // vals     : matrix (n_persons x n_items) with each row the valuation of that person for the items
  // alloc    : index of the person to which each item belongs
  // maxiter  : maximum number of iterations
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
  arma::mat avgval = get_avgval(vals, n_persons); // initialize average valuations
  arma::mat valmat = get_valmat(vals, alloc, n_items, n_persons); // get value of each set of items for each person
  arma::vec socvec = arma::zeros(maxiter + 1);    // initialize social inequality through iterations
  double minsoc = get_fnV(valmat, n_persons, avgval); // social inequality at initial stage
  socvec(0) = minsoc;
  
  // iterate
  int iter = 0;
  bool converged = false;
  int status = 1;
  while (iter < maxiter && !converged) {
    
    // sample items to give to a different owner
    int item = arma::randi(1, arma::distr_param(0, n_items - 1))(0);
    int addperson = arma::randi(1, arma::distr_param(1, n_persons - 1))(0);
    
    // try the reassignment and update if improvements are made
    arma::ivec alloctemp = alloc - 1;
    arma::mat valmattemp = valmat;
    int oldperson = alloctemp(item);
    int newperson = (oldperson + addperson) % n_persons;
    valmattemp.col(oldperson) -= vals.col(item);
    valmattemp.col(newperson) += vals.col(item);
    alloctemp(item) = newperson;
    double soctemp = get_fnV(valmattemp, n_persons, avgval);
    socvec(1 + iter) = soctemp;
    
    // update if lower maxenvy
    if (soctemp < minsoc) {
      minsoc = soctemp;
      alloc = alloctemp + 1;
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
