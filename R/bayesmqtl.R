# Wrapper function which both carries out mixture model estimation and variational inference.
bayesmqtl <- function(Y, X, list_hyper, list_init,
                      tol_mix = 0.1, tol_vb = 0.1,
                      maxit = 1000, verbose = FALSE) {

  # Control over numerical under/overflow
  eps <- 1e-6
  Y <- pmin(pmax(Y, eps), 1 - eps)

  # Run mixture modelling
  mix_model_fit <- fit_mixture_model_(Y, tol = tol_mix, obj_param = "all", sens_param = "none", true_params = NULL, sens_param_val = NULL, digamma_approx = FALSE) # mix_model_fit is a list where each element contains pi, alphas, betas, and the current log-likelihood

  finalit <- length(mix_model_fit)

  mix_model_fit <- mix_model_fit[[finalit]]
  filt_out <- mix_model_fit$filt_out

  mix_model_fit <- mix_model_fit$param_estimates[, -6]

  if (length(filt_out) > 0) {

    Y <- Y[, -filt_out]
    list_init$mu_gam_0_vb <- list_init$mu_gam_0_vb[-filt_out]
    list_init$mu_gam_1_vb <- list_init$mu_gam_1_vb[, -filt_out]
    list_init$tau_inv2_vb <- list_init$tau_inv2_vb[-filt_out]
  }
  else {

    filt_out <- "none"
  }

  # Prepare input for bayesmqtl_core_()
  d <- ncol(Y)

  pi_t <- mix_model_fit[, 1]
  alpha_0_t <- mix_model_fit[, 2]
  beta_0_t <- mix_model_fit[, 3]
  alpha_1_t <- mix_model_fit[, 4]
  beta_1_t <- mix_model_fit[, 5]

  logP_lower <- mapply(FUN = dbeta, x = as.data.frame(Y), shape1 = alpha_0_t, shape2 = beta_0_t, log = rep(TRUE, d))
  logP_upper <- mapply(FUN = dbeta, x = as.data.frame(Y), shape1 = alpha_1_t, shape2 = beta_1_t, log = rep(TRUE, d))

  # Run variational inference
  bayesmqtl <- bayesmqtl_core_(Y, X, logP_upper, logP_lower, list_hyper, list_init, tol = tol_vb, maxit, full_output = TRUE, verbose)

  bayesmqtl <- append(bayesmqtl, filt_out)
  names(bayesmqtl)[length(names(bayesmqtl))] <- "filt_out" # Need to include dplyr as a package dependency

  return(bayesmqtl)
}
