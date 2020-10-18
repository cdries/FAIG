#' Allocate indivisible goods
#'
#' wrapper function to allocate indivisible goods
#'
#'
#' There are currently four algorithms implemented: 1. randselect 2. localtrades 3. mincov and
#' 4. mincovtarget. For algorithms 1 and 2, there is the choice to use as objective function
#' either maxenvy (to be minimized) or social inequality (to be minimized). Algorithms 3 and 4
#' use the social inequality algorithm.
#' Control parameters include 'maxiter' for maximum number of iterations (default 1e5) and
#' 'eps' the tolerance to stop when V < eps (default 1e-6).
#'
#' @name allocate
#' @encoding UTF-8
#' @concept allocate
#' @param vals valuation matrix, each row represents the value for this agent for each of the items (columns)
#' @param algo algorithm, one of (mincov, mincovtarget, localtrades, randselect)
#' @param obj objective value to minimize, one of (soc, maxenvy, maxutility); only relevant for algorithms 
#' localtrades and randselect
#' @param alloc0 initial allocation, either a vector of length n_items containing the
#' index of the person to which each item belongs, or 'random', in which case we generate a random
#' initial allocation; not relevant when algo='randselect"
#' @param maxiter maximum number of iterations, default 1e5
#' @param maxnoimprove convergence criterium in number of steps yielding no improvement, default 1e3
#' @param eps absolute convergence criterion, default 1e-6
#' @param target target value or vector (length n_persons), only relevant for mincovtarget
#' @author Dries Cornilly
#' @references
#' Cornilly, D., Puccetti, G., Rüschendorf, L., & Vanduffel, S. (2020). 
#' New algorithms for fair allocation of indivisible goods with minimum inequality or minimum envy
#' criteria.
#'
#' @import Rcpp
#' @useDynLib FAIG
#' @export allocate
allocate <- function(vals, algo='mincov', obj='soc', alloc0='random', 
                     maxiter=1e5, maxnoimprove=1e3, eps=1e-6, target=0) {
  
  # initialize properties
  n_items  <- ncol(vals)
  n_persons <- nrow(vals)
  
  # call the requested algorithm
  if (algo == 'mincov') {
    out <- mincov_wrapper(vals, alloc0, n_items, n_persons, maxiter, maxnoimprove, eps)
  } else if (algo == 'mincovtarget') {
    out <- mincovtarget_wrapper(vals, alloc0, n_items, n_persons, maxiter, maxnoimprove, eps, target)
  } else if (algo == 'randselect') {
    out <- randselect_wrapper(vals, obj, maxiter, maxnoimprove, eps)
  } else if (algo == 'localtrades') {
    out <- localtrades_wrapper(vals, alloc0, obj, n_persons, n_items, maxiter, maxnoimprove, eps)
  } else {
    warning('Chosen algorithm not implemented.')
  }

  return (out) 
}
