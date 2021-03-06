# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @export get_valmat
get_valmat <- function(vals, alloc, n_items, n_persons) {
    .Call('_FAIG_get_valmat', PACKAGE = 'FAIG', vals, alloc, n_items, n_persons)
}

#' @export
get_maxenvy <- function(valmat, n) {
    .Call('_FAIG_get_maxenvy', PACKAGE = 'FAIG', valmat, n)
}

#' @export
get_avgval <- function(valmat, n) {
    .Call('_FAIG_get_avgval', PACKAGE = 'FAIG', valmat, n)
}

#' @export
get_fnV <- function(valmat, n, avgval) {
    .Call('_FAIG_get_fnV', PACKAGE = 'FAIG', valmat, n, avgval)
}

#' @export
get_util <- function(valmat) {
    .Call('_FAIG_get_util', PACKAGE = 'FAIG', valmat)
}

testfunc <- function(oldperson, addperson, n_persons) {
    .Call('_FAIG_testfunc', PACKAGE = 'FAIG', oldperson, addperson, n_persons)
}

localtrades_envy <- function(vals, alloc, maxiter, maxnoimprove, eps) {
    .Call('_FAIG_localtrades_envy', PACKAGE = 'FAIG', vals, alloc, maxiter, maxnoimprove, eps)
}

localtrades_social <- function(vals, alloc, maxiter, maxnoimprove, eps) {
    .Call('_FAIG_localtrades_social', PACKAGE = 'FAIG', vals, alloc, maxiter, maxnoimprove, eps)
}

localtrades_utility <- function(vals, alloc, maxiter, maxnoimprove, eps) {
    .Call('_FAIG_localtrades_utility', PACKAGE = 'FAIG', vals, alloc, maxiter, maxnoimprove, eps)
}

mincov <- function(vals, alloc, beta, maxiter, maxnoimprove, eps) {
    .Call('_FAIG_mincov', PACKAGE = 'FAIG', vals, alloc, beta, maxiter, maxnoimprove, eps)
}

mincovtarget <- function(vals, alloc, beta, target, maxiter, maxnoimprove, eps) {
    .Call('_FAIG_mincovtarget', PACKAGE = 'FAIG', vals, alloc, beta, target, maxiter, maxnoimprove, eps)
}

random_alloc <- function(n_items, n_persons) {
    .Call('_FAIG_random_alloc', PACKAGE = 'FAIG', n_items, n_persons)
}

randselect_envy <- function(vals, maxiter, maxnoimprove, eps) {
    .Call('_FAIG_randselect_envy', PACKAGE = 'FAIG', vals, maxiter, maxnoimprove, eps)
}

randselect_social <- function(vals, maxiter, maxnoimprove, eps) {
    .Call('_FAIG_randselect_social', PACKAGE = 'FAIG', vals, maxiter, maxnoimprove, eps)
}

