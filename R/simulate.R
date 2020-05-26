#' Simulation scenario used in the paper
#'
#' simulate standardized valuation matrices
#'
#' @name simulate
#' @encoding UTF-8
#' @concept simulate
#' @param n number of agents
#' @param d number of items; should be a multiple of the number of agents
#' @param eps parameter governing the dependence between the agents' valuations
#' @author Dries Cornilly
#' @references
#' Cornilly, D., Puccetti, G., Rüschendorf, L., & Vanduffel, S. (2020). 
#' New algorithms for fair allocation of indivisible goods with minimum inequality or minimum envy
#' criteria.
#'
#' @importFrom stats runif
#' @export simulate
simulate <- function(n, d, eps=0.5) {
  
  # generate values first agent and ranges for the other agents
  x <- matrix(stats::runif(d * n) * 100, nrow = n, ncol = d)
  mm <- x[1,] * (1 - eps)
  MM <- x[1,] * (1 + eps)
  
  # generate values for the other agents
  for (ii in 2:n) {
    x[ii,] <- runif(d, min = mm, max = MM)
  }
  
  # make no-envy solution - iterate over the agents
  ipp <- round(d / n)
  for (ii in 1:n) {
    values <- sapply(1:n, function(a) sum(x[ii, ((a - 1) * ipp + 1):(a * ipp)]))
    indmax <- which.max(values)
    if (indmax != ii) {
      temp <- x[ii, ((ii - 1) * ipp + 1):(ii * ipp)]
      x[ii, ((ii - 1) * ipp + 1):(ii * ipp)] <- x[ii, ((indmax - 1) * ipp + 1):(indmax * ipp)]
      x[ii, ((indmax - 1) * ipp + 1):(indmax * ipp)] <- temp
    }
  }
  
  # normalize
  x <- x / rowSums(x) * 100
  
  # randomize
  vals <- x[, sample(1:d, d, replace = FALSE)]
  
  return (vals)
}

