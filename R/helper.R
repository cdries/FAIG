
#' @export ppl_1n
ppl_1n <- function(valmat) {
  k <- nrow(valmat)
  tmp <- sum(diag(valmat) / rowSums(valmat) > 1 / k - 1e-10)
  return (tmp)
}


#' @export ppl_max
ppl_max <- function(valmat) {
  tmp <- sum(diag(valmat) - apply(valmat, 1, max) > -1e-10)
  return (tmp)
}


get_beta <- function(vals, n_persons, n_items) {
  
  # initialize
  beta <- matrix(1, nrow = n_persons, ncol = n_items)
  
  # iterate over persons
  for (ii in 2:n_persons) {
    beta[ii,] <- vals[ii,] / vals[1,]
  }
  
  # set non-finite values (due to division by zero) to zero
  beta[!is.finite(beta)] <- 0
  
  return (beta)
}


mincov_wrapper <- function(vals, alloc0, n_items, n_persons, maxiter, eps) {
  
  # get initial allocation if necessary
  if (alloc0[1] == 'random') {
    alloc0 <- c(random_alloc(n_items, n_persons))
  }
  
  # get beta
  beta <- get_beta(vals, n_persons, n_items)
  
  # call mincov
  out <- mincov(vals, alloc0, beta, maxiter, eps)
  
  return (out)
}


mincovtarget_wrapper <- function(vals, alloc0, n_items, n_persons, maxiter, eps, target) {
  
  # get initial allocation if necessary
  if (alloc0[1] == 'random') {
    alloc0 <- c(random_alloc(n_items, n_persons))
  }
  
  # get beta
  beta <- get_beta(vals, n_persons, n_items)
  
  # call mincov
  if (length(target) == 1) {
    target <- rep(target, n_persons)
  }
  out <- mincovtarget(vals, alloc0, beta, target, maxiter, eps)
  
  return (out)
}


randselect_wrapper <- function(vals, obj, maxiter, eps) {
  
  # call randselect implementation depending on the objective
  if (obj == 'soc') {
    out <- randselect_social(vals, maxiter, eps)
  } else if (obj == 'maxenvy') {
    out <- randselect_envy(vals, maxiter, eps)
  } else {
    warning('Objective not implemented.')
  }
  
  return (out)
}


localtrades_wrapper <- function(vals, alloc0, obj, n_persons, n_items, maxiter, eps) {
  
  # get initial allocation if necessary
  if (alloc0[1] == 'random') {
    alloc0 <- c(random_alloc(n_items, n_persons))
  }
  
  # call localtrades implementation depending on the objective
  if (obj == 'soc') {
    out <- localtrades_social(vals, alloc0, maxiter, eps)
  } else if (obj == 'maxenvy') {
    out <- localtrades_envy(vals, alloc0, maxiter, eps)
  } else {
    warning('Objective not implemented')
  }
  
  return (out)
}
