# Wrapper function which both carries out mixture model estimation and variational inference.
bayesmqtl <- function(Y, X, list_hyper, list_init,
                      tol_mix = 0.1, tol_vb = 0.1,
                      maxit = 1000, verbose = FALSE) {

  d <- ncol(Y)

  # Run mixture modelling
  mix_model_fit <- fit_mixture_model_(Y, tol = tol_mix, obj_param = "all", sens_param = "none", true_params = NULL, digamma_approx = TRUE) # mix_model_fit is a list where each element contains pi, alphas, betas, and the current log-likelihood

  finalit <- length(mix_model_fit)
  mix_model_fit <- mix_model_fit[[finalit]]
  mix_model_fit <- mix_model_fit$param_estimates[, -6]

  # Prepare input for bayesmqtl_core_()
  pi_t <- mix_model_fit[, 1]
  alpha_0_t <- mix_model_fit[, 2]
  beta_0_t <- mix_model_fit[, 3]
  alpha_1_t <- mix_model_fit[, 4]
  beta_1_t <- mix_model_fit[, 5]

  logP_lower <- mapply(FUN = dbeta, x = as.data.frame(Y), shape1 = alpha_0_t, shape2 = beta_0_t, log = rep(TRUE, d))
  logP_upper <- mapply(FUN = dbeta, x = as.data.frame(Y), shape1 = alpha_1_t, shape2 = beta_1_t, log = rep(TRUE, d))

  # Run variational inference
  bayesmqtl <- bayesmqtl_core_(Y, X, logP_upper, logP_lower, list_hyper, list_init, tol = tol_vb, maxit, full_output = FALSE, verbose)

  return(bayesmqtl)
}
