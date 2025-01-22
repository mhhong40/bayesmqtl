#' Internal core function to call the variational algorithm.
#' Structure is largely adapted from Hélène's locus core code.

bayesmqtl_core_ <- function(z, X, lambda, eta, kappa, tau,
                            theta, sig2_theta,
                            tol = 0.1, maxit = 1000, anneal = NULL,
                            checkpoint_path = NULL, verbose = TRUE) {

  n <- nrow(z)
  d <- ncol(z)
  p <- ncol(X)

  eps <- .Machine$double.eps^0.5

  # Preparing annealing, if any...
  if (is.null(anneal)) {

    annealing <- FALSE
    c <- 1
  }
  else {

    annealing <- TRUE
    ladder <- get_annealing_ladder_(anneal, verbose)
    c <- ladder[1]
  }

  tau_vb <- tau
  theta_vb <- mu_theta_vb <- theta # two objects are created to make the computations in C++ Eigen easier
  sig2_theta_vb <- matrix(rep(sig2_theta, d), ncol = d)

  m2_theta <- update_m2_theta_vb_(mu_theta_vb, sig2_theta_vb)

  mat_x_m1 <- X %*% theta_vb # mat_x_m1 contains the predicted values for each z
                             # It's a matrix with dimensions n x d.

  converged <- FALSE
  lb_new <- -Inf
  it <- 0

  # Iterative algorithm
  while((!converged) & (it < maxit)) {

    lb_old <- lb_new
    it <- it + 1

    # % #
    eta_vb <- update_eta_vb_(n, eta, c = c)
    kappa_vb <- update_kappa_vb_(z, kappa, mat_x_m1, m2_theta, theta_vb, c = c)

    tau_vb <- eta_vb / kappa_vb # this somehow became negative
    # % #

    # % #
    lambda_vb <- update_lambda_vb_(lambda, m2_theta, c)
    sig2_vb <- matrix(rep(1 / lambda_vb, d), ncol = d)
    # % #

    # % #
    sig2_theta_vb <- update_sig2_theta_vb_(n, tau_vb, sig2_vb, c)

    for (s in 1:p) {

      mat_x_m1 <- mat_x_m1 - tcrossprod(X[, s], theta_vb[s, ]) # take away the contribution of the sth predictor

      mu_theta_vb[s, ] <- c * sig2_theta_vb[s, ] * tau_vb *
        (t(z - mat_x_m1) %*% X[, s])

      theta_vb[s, ] <- mu_theta_vb[s, ]

      mat_x_m1 <- mat_x_m1 + tcrossprod(X[, s], theta_vb[s, ])
    }

    ## TO DO: add performance optimizer flags to perform computationally expensive subroutines
    ##        in C++ using the coreLoop() function below:
    # coreLoop(z, X, mat_x_m1, theta_vb,
    #          mu_theta_vb, tau_vb, sig2_theta_vb, c = c) # Updates theta_vb
    # % #

    m2_theta <- update_m2_theta_vb_(mu_theta_vb, sig2_theta_vb)

    # % #
    # Adapted from Hélène's locus>locus_core code.
    if (annealing) {

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste0("Temperature = ", format(1 / c, digits = 4), "\n\n"))

      c <- ifelse(it < length(ladder), ladder[it + 1], 1)

      if (isTRUE(all.equal(c, 1))) {

        annealing <- FALSE

        if (verbose)
          cat("** Exiting annealing mode. **\n\n")
      }
    }
    else {

      log_sig2_theta_vb <- log(sig2_theta_vb + eps) - 1/2
      log_tau_vb <- log(tau_vb + eps)

      lb_new <- elbo_(n, z, mat_x_m1, theta_vb, kappa, kappa_vb, log_tau_vb, m2_theta, tau_vb,
                      eta, eta_vb, sig2_theta_vb, log_sig2_theta_vb, lambda, lambda_vb, sig2_vb)

      if (verbose)
        cat(paste0("ELBO = ", format(lb_new), "\n\n",
                   "Difference in ELBO = ", format(abs(lb_new - lb_old)), "\n\n"))

      # if (lb_new + eps < lb_old)
      #  stop("ELBO not increasing monotonically. Exit. ")

      converged <- (abs(lb_new - lb_old) < tol)

      # checkpoint_(it, checkpoint_path, theta_vb, converged, lb_new, lb_old)

    }

    checkpoint_clean_up_(checkpoint_path)

    if(verbose) {

      if(converged) {

        cat(paste0("Convergence obtained after ", format(it), " interations. \n",
                   "Optimal marginal log-likelihood variational lower bound ",
                   "(ELBO) = ", format(lb_new), ". \n\n"))
      }
      # else {

      #  warning("Maximal number of iterations reached before convergence. Exit.")
      # }
    }

  }

  # When convergence has been reached:
  lb_opt <- lb_new

  names_x <- colnames(X)
  names_z <- colnames(z)

  rownames(theta_vb) <- rownames(sig2_theta_vb) <- names_x
  colnames(theta_vb) <- colnames(sig2_theta_vb) <- names_z

  diff_lb <- abs(lb_opt - lb_old)

  annealing <- ifelse(is.null(anneal), FALSE, anneal[1])

  create_named_list_(theta_vb, sig2_theta_vb,
                     converged, it, lb_opt, diff_lb, annealing)

}

elbo_ <- function(n, z, mat_x_m1, theta_vb, kappa, kappa_vb, log_tau_vb, m2_theta, tau_vb,
                  eta, eta_vb, sig2_theta_vb, log_sig2_theta_vb, lambda, lambda_vb, sig2_vb) {

  eps <- .Machine$double.eps^0.5

  n <- nrow(z)

  eta_vb <- update_eta_vb_(n, eta)
  kappa_vb <- update_kappa_vb_(z, kappa, mat_x_m1, theta_vb, m2_theta)
  log_tau_vb <- digamma(eta_vb) - log(kappa_vb + eps) # <- this is what throws the error -- log(kappa_vb + eps) yields NaNs

  log_sig2_theta_vb <- log(sig2_theta_vb) - 1/2 # approximating E(log(sig2_theta))

  lambda_vb <- update_lambda_vb_(lambda, m2_theta)

  elbo_A <- e_z_(n, kappa, kappa_vb, log_tau_vb, tau_vb, m2_theta)

  elbo_B <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)

  elbo_C <- e_theta_(sig2_theta_vb, log_sig2_theta_vb, m2_theta)

  elbo_D <- e_sig2_(lambda, lambda_vb, sig2_vb)

  elbo_A + elbo_B + elbo_C + elbo_D

}
