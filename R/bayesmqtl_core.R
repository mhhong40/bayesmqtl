# Core function containing the iterative variational algorithm for bayesmqtl
bayesmqtl_core_ <- function(Y, X, logP_upper, logP_lower, list_hyper, list_init, tol, maxit, full_output = TRUE, verbose) { # logP_upper and logP_lower are computed outside of the core function and supplied to it

  eps <- .Machine$double.eps^0.5

  n <- nrow(Y)
  p <- ncol(X)
  d <- ncol(Y)

  # Extract initialized model parameters
  mu_gam_0_vb <- list_init$mu_gam_0_vb
  mu_gam_1_vb <- list_init$mu_gam_1_vb
  tau_inv2_vb <- list_init$tau_inv2_vb
  lambda_inv2_vb <- list_init$lambda_inv2_vb

  # Extract model hyperparameters
  eta_0 <- list_hyper$eta_0
  sig_0_inv2 <- list_hyper$sig_0_inv2

  converged <- FALSE
  lb_new_vec <- rep(-Inf, 8) # 8 model parameter types
  lb_new <- -Inf
  elbos <- c()
  it <- 0

  sig2_gam_0_vb <- update_sig2_gam_0_(sig_0_inv2, n, d) # Constant, so can be placed outside of the loop
  sig2_gam_1_vb <- update_sig2_gam_1_(X, tau_inv2_vb, lambda_inv2_vb)

  xi_vb <- mu_gam_0_vb + X %*% mu_gam_1_vb

  log_Phi_xi_vb <- pnorm(xi_vb, log.p = TRUE)
  log_1_Phi_xi_vb <- pnorm(xi_vb, log.p = TRUE, lower.tail = FALSE)

  z_vb <- update_z_(logP_upper, logP_lower, log_Phi_xi_vb, log_1_Phi_xi_vb)
  rho_vb <- update_rho_(X, z_vb, xi_vb)

  # Preset object
  resid_mu_gam_1 <- matrix(nrow = n, ncol = d)

  # Basic code profiling
  start_time <- proc.time()[[3]]

  # Track individual changes in ELBO
  lb_old_vec <- c()
  lb_diffs <- as.matrix(t(lb_new_vec))
  colnames(lb_diffs) <- c("y", "z, rho", "gam_0", "gam_1", "lambda", "tau", "a", "b")

  while((!converged) & (it < maxit)) {

    lb_old_vec <- lb_new_vec
    lb_old <- lb_new
    it <- it + 1

    if (verbose)
      cat(paste0("Iteration ", format(it), "... \n"))

    # % #
    a_inv_vb <- update_a_inv_(lambda_inv2_vb, p, log = FALSE)
    b_inv_vb <- update_b_inv_(tau_inv2_vb, d, log = FALSE)

    eta_lambda <- update_eta_lambda_(tau_inv2_vb, sig2_gam_1_vb, mu_gam_1_vb, a_inv_vb, d)
    eta_tau <- update_eta_tau_(lambda_inv2_vb, sig2_gam_1_vb, mu_gam_1_vb, b_inv_vb, p)

    lambda_inv2_vb <- update_lambda_inv2_(eta_lambda, p, d, log = FALSE)
    tau_inv2_vb <- update_tau_inv2_(eta_tau, p, d, log = FALSE)
    # % #

    # % #
    sig2_gam_0_vb <- update_sig2_gam_0_(sig_0_inv2, n, d)
    sig2_gam_1_vb <- update_sig2_gam_1_(X, tau_inv2_vb, lambda_inv2_vb)

    mu_gam_0_vb <- update_mu_gam_0_(sig_0_inv2, sig2_gam_0_vb, eta_0, X, mu_gam_1_vb, rho_vb)
    mu_gam_1_vb <- update_mu_gam_1_(sig2_gam_1_vb, X, d, rho_vb, mu_gam_0_vb, mu_gam_1_vb, resid_mu_gam_1, vectorized = TRUE)
    # % #

    # % #
    xi_vb <- mu_gam_0_vb + X %*% mu_gam_1_vb # xi_vb gets updated

    log_Phi_xi_vb <- pnorm(xi_vb, log.p = TRUE)
    log_1_Phi_xi_vb <- pnorm(xi_vb, log.p = TRUE, lower.tail = FALSE)

    z_vb <- update_z_(logP_upper, logP_lower, log_Phi_xi_vb, log_1_Phi_xi_vb)
    rho_vb <- update_rho_(X, z_vb, xi_vb)
    # % #

    lb_new_vec <- elbo_(logP_upper, logP_lower, Y, X, z_vb, log_Phi_xi_vb, log_1_Phi_xi_vb,
                    sig2_gam_0_vb, sig2_gam_1_vb,
                    sig_0_inv2, eta_0, mu_gam_0_vb, mu_gam_1_vb,
                    tau_inv2_vb, lambda_inv2_vb,
                    eta_tau, eta_lambda, a_inv_vb, b_inv_vb)

    lb_diff_vec <- lb_new_vec - lb_old_vec
    lb_diffs <- rbind(lb_diffs, lb_diff_vec)

    lb_new <- sum(lb_new_vec)
    elbos <- append(elbos, lb_new)

    cat("ELBO total: ", lb_new, "\n")

    if (lb_new + eps < lb_old)
      stop("ELBO not increasing monotonically. Exit.")

    if (lb_new - lb_old + eps < tol)
      converged <- TRUE
  }

  if (verbose != 0) {
    if (converged) {
      cat(paste0("Convergence obtained after ", format(it), " iterations. \n",
                 "Optimal marginal log-likelihood variational lower bound ",
                 "(ELBO) = ", format(lb_new), ". \n\n"))
    } else {
      warning("Maximal number of iterations reached before convergence. Exit.")
    }
  }
  lb_opt <- lb_new

  diff_lb <- lb_opt - lb_old

  end_time <- proc.time()[[3]]
  time_elapsed <- end_time - start_time

  if(full_output) { # Only for internal use

    output <- create_named_list_(z_vb, rho_vb, sig2_gam_0_vb, sig2_gam_1_vb,
                                 mu_gam_0_vb, mu_gam_1_vb, tau_inv2_vb, lambda_inv2_vb,
                                 a_inv_vb, b_inv_vb, elbos, lb_diffs, time_elapsed)

    return(output)
  }
  else {

    names_x <- colnames(X)
    names_y <- colnames(Y)
    names_n <- rownames(Y)

    rownames_z_vb <- rownames_rho_vb <- names_n
    colnames_z_vb <- colnames_rho_vb <- names_y

    names(mu_gam_0_vb) <- names(mu_gam_1_vb) <- names_y

    create_named_list_(z_vb, rho_vb, mu_gam_0_vb, mu_gam_1_vb,
                       n, p, d, converged, it, maxit, tol, lb_opt, diff_lb, lb_diffs, time_elapsed)
  }
}

# ELBO function (wrapper for all the ELBO terms coded in elbo.R)
elbo_ <- function(logP_upper, logP_lower, Y, X, z_vb, log_Phi_xi_vb, log_1_Phi_xi_vb,
                  sig2_gam_0_vb, sig2_gam_1_vb,
                  sig_0_inv2, eta_0, mu_gam_0_vb, mu_gam_1_vb,
                  tau_inv2_vb, lambda_inv2_vb,
                  eta_tau, eta_lambda, a_inv_vb, b_inv_vb) {

  n <- nrow(Y)
  d <- ncol(Y)
  p <- ncol(X)

  # Make sure params are up-to-date with log-expectations
#  sig2_gam_0_vb <- update_sig2_gam_0_(sig_0_inv2, n, d)
#  sig2_gam_1_vb <- update_sig2_gam_1_(X, tau_inv2_vb, lambda_inv2_vb)

#  eta_lambda <- update_eta_lambda_(tau_inv2_vb, sig2_gam_1_vb, mu_gam_1_vb, a_inv_vb, d)
#  eta_tau <- update_eta_tau_(lambda_inv2_vb, sig2_gam_1_vb, mu_gam_1_vb, b_inv_vb, p)

  # Expectations of logs only needed for computing the ELBO
  log_a_inv_vb <- update_a_inv_(lambda_inv2_vb, p, log = TRUE)
  log_b_inv_vb <- update_b_inv_(tau_inv2_vb, d, log = TRUE)

  log_tau_inv2_vb <- update_tau_inv2_(eta_tau, p, d, log = TRUE)
  log_lambda_inv2_vb <- update_lambda_inv2_(eta_lambda, p, d, log = TRUE)

  d_a <- lambda_inv2_vb + 1
  d_b <- tau_inv2_vb + 1

  # Actual ELBO terms
  elbo_A <- elbo_y_(logP_upper, logP_lower, z_vb)

  elbo_B <- elbo_z_rho_(X, z_vb, log_Phi_xi_vb, log_1_Phi_xi_vb, sig2_gam_0_vb, sig2_gam_1_vb)

  elbo_C <- elbo_gam_0_(sig_0_inv2, eta_0, mu_gam_0_vb, sig2_gam_0_vb)

  elbo_D <- elbo_gam_1_(log_tau_inv2_vb, log_lambda_inv2_vb, tau_inv2_vb, lambda_inv2_vb, mu_gam_1_vb, sig2_gam_1_vb)

  elbo_E <- elbo_lambda_(a_inv_vb, log_a_inv_vb, log_lambda_inv2_vb, eta_lambda, lambda_inv2_vb, d)

  elbo_F <- elbo_tau_(b_inv_vb, log_b_inv_vb, log_tau_inv2_vb, eta_tau, tau_inv2_vb, p)

  elbo_G <- elbo_a_(d_a, a_inv_vb, log_a_inv_vb)

  elbo_H <- elbo_b_(d_b, b_inv_vb, log_b_inv_vb)

#  cat("elbo_y:", elbo_A, "\n")
#  cat("elbo_z_rho:", elbo_B, "\n")
#  cat("elbo_gam_0:", elbo_C, "\n")
#  cat("elbo_gam_1:", elbo_D, "\n")
#  cat("elbo_lambda:", elbo_E, "\n")
#  cat("elbo_tau:", elbo_F, "\n")
#  cat("elbo_a:", elbo_G, "\n")
#  cat("elbo_b:", elbo_H, "\n")

  return(c(elbo_A, elbo_B, elbo_C, elbo_D, elbo_E, elbo_F, elbo_G, elbo_H))
}
