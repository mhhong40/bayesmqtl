# Functions comprising the ELBO terms in elbo() from bayesmqtl_core.R

elbo_y_ <- function(logP_upper, logP_lower, z_vb) {

  elbo_y <- z_vb * logP_upper + (1 - z_vb) * logP_lower
  elbo_y <- sum(elbo_y)

  return(elbo_y)
}

elbo_z_rho_ <- function(X, z_vb, log_Phi_xi_vb, log_1_Phi_xi_vb, sig2_gam_0_vb, sig2_gam_1_vb) {

  elbo_no_entropy_z <- sweep(z_vb * log_Phi_xi_vb + (1 - z_vb) * log_1_Phi_xi_vb - 1/2 * (X^2 %*% sig2_gam_1_vb), 2, 1/2 * sig2_gam_0_vb, "-")

  entropy_z <- - (z_vb * log(z_vb) + (1 - z_vb) * log(1 - z_vb))
  entropy_z[is.nan(entropy_z)] <- 0

  elbo_z_rho <- sum(elbo_no_entropy_z) + sum(entropy_z)

  return(elbo_z_rho)
}

elbo_gam_0_ <- function(sig_0_inv2, eta_0, mu_gam_0_vb, sig2_gam_0_vb) {

  elbo_gam_0 <- 1/2 * (-log(sig_0_inv2) - sig_0_inv2 * ((mu_gam_0_vb - eta_0)^2 + sig2_gam_0_vb) + log(sig2_gam_0_vb) + 1)
  elbo_gam_0 <- sum(elbo_gam_0)

  return(elbo_gam_0)
}

elbo_gam_1_ <- function(log_tau_inv2_vb, log_lambda_inv2_vb, tau_inv2_vb, lambda_inv2_vb, mu_gam_1_vb, sig2_gam_1_vb) {

  elbo_gam_1 <- 1/2 * (-outer(log_lambda_inv2_vb, log_tau_inv2_vb, "+") + outer(lambda_inv2_vb, tau_inv2_vb, "*") * (mu_gam_1_vb^2 + sig2_gam_1_vb) + log(sig2_gam_1_vb) + 1) # Functional version
  # elbo_gam_1 <- 1/2 * (outer(log_lambda_inv2_vb, log_tau_inv2_vb, "+") - outer(lambda_inv2_vb, tau_inv2_vb, "*") * (mu_gam_1_vb^2 + sig2_gam_1_vb) + log(sig2_gam_1_vb) + 1) # Derived expression

  elbo_gam_1 <- sum(elbo_gam_1)

  return(elbo_gam_1)
}

elbo_lambda_ <- function(a_inv_vb, log_a_inv_vb, log_lambda_inv2_vb, eta_lambda, lambda_inv2_vb, d) {

  elbo_lambda <- 1/2 * ((eta_lambda - a_inv_vb)*lambda_inv2_vb - d*log_lambda_inv2_vb + log_a_inv_vb - (d+1)*log(eta_lambda)) # Functional version
  # elbo_lambda <- 1/2 * (d*log_lambda_inv2_vb + log_a_inv_vb - (d+1)*log(eta_lambda)) + (eta_lambda - a_inv_vb)*lambda_inv2_vb # Derived expression

  elbo_lambda <- sum(elbo_lambda)

  return(elbo_lambda)
}

elbo_tau_ <- function(b_inv_vb, log_b_inv_vb, log_tau_inv2_vb, eta_tau, tau_inv2_vb, p) {

  elbo_tau <- 1/2 * ((eta_tau - b_inv_vb)*tau_inv2_vb - p*log_tau_inv2_vb + log_b_inv_vb - (p+1)*log(eta_tau)) # Functional version
  # elbo_tau <- 1/2 * (p*log_tau_inv2_vb + log_b_inv_vb - (p+1)*log(eta_tau)) + (eta_tau - b_inv_vb)*tau_inv2_vb # Derived expression

  elbo_tau <- sum(elbo_tau)

  return(elbo_tau)
}

elbo_a_ <- function(d_a, a_inv_vb, log_a_inv_vb) {

  elbo_a <- (d_a - 1) * a_inv_vb - 1/2 * log_a_inv_vb - log(d_a)
  elbo_a <- sum(elbo_a)

  return(elbo_a)
}

elbo_b_ <- function(d_b, b_inv_vb, log_b_inv_vb) {

  elbo_b <- (d_b - 1) * b_inv_vb - 1/2 * log_b_inv_vb - log(d_b)
  elbo_b <- sum(elbo_b)

  return(elbo_b)
}
