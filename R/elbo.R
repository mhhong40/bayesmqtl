# Functions comprising the ELBO terms in elbo() from bayesmqtl_core.R

elbo_y_ <- function(logP, log_1P, z_vb) {

  elbo_y <- z_vb * logP + (1 - z_vb) * log_1P
  elbo_y <- sum(elbo_y)

  return(elbo_y)
}

elbo_z_rho_ <- function(X, z_vb, log_Phi_xi_vb, log_1_Phi_xi_vb, sig2_gam_0_vb, sig2_gam_1_vb) {

  elbo_z_rho <- z_vb * log_Phi_xi_vb + (1 - z_vb) * log_1_Phi_xi_vb - 1/2 * sig2_gam_0_vb - 1/2 * sweep(X^2, 2, sig2_gam_1_vb, "*")
  elbo_z_rho <- sum(elbo_z_rho)

  return(elbo_z_rho)
}

elbo_gam_0_ <- function(sig_0_inv2, eta_0, mu_gam_0_vb, sig2_gam_0_vb) {

  elbo_gam_0 <- 1/2 * (-log(1/sig_0_inv2) - sig_0_inv2 * (mu_gam_0_vb - eta_0)^2 - sig_0_inv2 * sig2_gam_0_vb + log(sig2_gam_0_vb))
  elbo_gam_0 <- sum(elbo_gam_0)

  return(elbo_gam_0_)
}

elbo_gam_1_ <- function(log_tau_inv2_vb, log_lambda_inv2_vb, tau_inv2_vb, lambda_inv2_vb, mu_gam_1_vb, sig2_gam_1_vb) {

  elbo_gam_1 <- 1/2 * (log_tau_inv2_vb + log_lambda_inv2_vb - tau_inv2_vb * lambda_inv2_vb * (mu_gam_1_vb^2 + sig2_gam_1_vb) + log(sig2_gam_1_vb))
  elbo_gam_1 <- sum(elbo_gam_1)

  return(elbo_gam_1)
}

elbo_lambda_ <- function(log_b_vb, log_lambda_inv2_vb, eta_lambda, lambda_inv2_vb) {

  elbo_lambda <- 1/2 * (log_b_vb - log_lambda_inv2_vb - log(eta_lambda) + log(lambda_inv2_vb))
  elbo_lambda <- sum(elbo_lambda)

  return(elbo_lambda)
}

elbo_tau_ <- function(log_a_vb, log_tau_inv2_vb, eta_tau, tau_inv2_vb) {

  elbo_tau <- 1/2 * (log_a_vb - log_tau_inv2_vb - log(eta_tau) + log(tau_inv2_vb))

  return(elbo_tau)
}

elbo_a_ <- function(d_a, a_vb, log_a_vb) {

  elbo_a <- (d_a - 1) * a_vb - 1/2 * log_a_vb - log(d_a)

  return(elbo_a)
}

elbo_b_ <- function(d_b, b_vb, log_b_vb) {

  elbo_b <- (d_b - 1) * b_vb - 1/2 * log_b_vb - log(d_b)
  elbo_b <- sum(elbo_b)

  return(elbo_b)
}
