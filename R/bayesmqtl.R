#' Contains code for the bayesmqtl() wrapper that executes
#' the full variational inference procedure.
#' Much of it is adapted from Hélène's locus>locus locus() function.
#' TO DO: write a proper description.
#' @export
bayesmqtl <- function(Y, X, z, bin_ind, mix_prop, alpha_0, beta_0, alpha_1, beta_1, list_hyper, list_init,
                      user_seed = NULL, tol = 0.1, maxit = 1000,
                      anneal = NULL, save_hyper = FALSE, save_init = FALSE,
                      full_output = FALSE, verbose = TRUE, checkpoint_path = NULL) {

  check_structure_(verbose, "vector", "logical", 1)

  # % #
  if (verbose) cat("== Preparing the data ... \n\n")

  check_annealing_(anneal)

  dat <- prepare_data_(Y, X, user_seed, tol, maxit, verbose, checkpoint_path)

  X <- dat$X
  Y <- dat$Y

  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(Y)
  k <- length(bin_ind)

  names_x <- colnames(X)
  names_y <- colnames(Y)

  if (verbose) cat("... done. == \n\n")
  # % #


  # % #
  if (verbose) cat("== Preparing the fixed shape parameters ... \n\n")

  bin_ind <- prepare_bin_ind_(bin_ind, Y, verbose)
  list_shape <- prepare_list_shape_(list_shape, Y, k, verbose)

  if (verbose) cat("... done. == \n\n")
  # % #


  # % #
  if (verbose) cat("== Preparing the hyperparameters ... \n\n")

  list_hyper <- prepare_list_hyper_(list_hyper, Y, p, names_x, names_y, verbose)

  if (verbose) cat("... done. == \n\n")

  if (verbose) cat("== Preparing the parameter initialization ... \n\n")

  list_init <- prepare_list_init_(list_init, Y, p, user_seed, verbose)

  if (verbose) cat("... done. == \n\n")
  # % #

  # % #
  if (verbose){
    cat(paste0("=================================================================== \n",
               "== Variational inference for class membership linear regressions == \n",
               "=================================================================== \n\n"))
  }

  sig2_theta <- 1/list_hyper$lambda

  # Initial hyper/parameter values taken out from list_hyper/init and passed into core.
  vb <- bayesmqtl_core_(z, X, list_hyper$lambda, list_hyper$eta, list_hyper$kappa,
                        list_init$tau, theta, sig2_theta,
                        tol = 0.1, maxit = 1000, anneal, verbose = TRUE)

  # % #

  # % #
  if(save_hyper) vb$list_hyper <- list_hyper
  if(save_init) vb$list_init <- list_init

  class(vb) <- "bayesmqtl"

  if(verbose) cat("... done. == \n\n")
  # % #

  vb
}
