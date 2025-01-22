# This file is part of the `bayesmqtl` package:
#   https://github.com/mhhong40/bayesmqtl
#
#' Gather model hyperparameters provided by the user.
#'
#' Currently, this function must be used to provide hyperparameter values for the model
#' used in \code{\link{bayesmqtl}}. A default setting has not been implemented yet,
#' but is on the to-do list.
#'
#' @param d Number of CpGs.
#' @param p Number of SNPs.
#' @param X A matrix of size n x p containing the SNP genotype data. Values must represent the minor allele dosage at each locus.
#'    Calculation of each SNP's mean and variance assumes that the SNPs follow Hardy-Weinberg equilibrium,
#'    i.e., are distributed along \eqn{Binom(2, \text{maf})}, where \eqn{\text{maf}} = the minor allele frequency.
#' @return An object of class "\code{hyper}" preparing user hyperparameters in a
#'    form passable to the \code{\link{bayesmqtl}} function.
#' @export
set_hyper <- function(d, p, X, kappa, eta, lambda) {

  # TO DO: code check_structure_ in utils.R
  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)
  d_hyper <- d

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)
  p_hyper <- p

  check_structure_(kappa, "vector", "double", c(1, d))
  if (length(kappa) == 1) kappa <- rep(kappa, d)

  check_structure_(eta, "vector", "double", c(1, d))
  if (length(eta) == 1) eta <- rep(eta, d)

  check_structure_(lambda, "vector", "double", c(1, p))
  if (length(lambda) == 1) lambda <- rep(lambda, p)

  list_hyper <- c(create_named_list_(d_hyper, p_hyper,
                                     kappa, eta, lambda))

  class(list_hyper) <- "hyper"

  list_hyper
}

#' Gather initial variational parameters provided by the user.
#'
#' Currently, this function must be used to provide initial values for the variational
#' parameters used in \code{\link{bayesmqtl}}. A default setting has not been implemented yet,
#' but is on the to-do list.
#'
#' @param d Number of CpGs.
#' @param p Number of SNPs.
#' @param theta Matrix of size p x d with initial values for the slope coefficients associated with each
#'    predictor SNP in the latent class membership. These values are updated with each iteration of the VI algorithm.
#' @param tau Vector of length d containing error variances
#' @return An object of class "\code{init}" preparing user initial values for
#'   the variational parameters in a form passable to the
#'   \code{\link{bayesmqtl}} function.
#' @export
set_init <- function(d, p, theta, tau) {

  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)
  d_init <- d

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)
  p_init <- p

  check_structure_(theta, "matrix", "double", c(p, d))

  check_structure_(tau, "vector", "double", c(1, d))
  if (length(tau) == 1) tau <- rep(tau, d)

  list_init <- create_named_list_(d_init, p_init, theta, tau)

  class(list_init) <- "init"

  list_init
}
