# Currently implemented in the multiple regression setting (p = 1 SNP, all d CpGs)

update_z_ <- function(logP, log_1P, log_Phi_xi_vb, log_1_Phi_xi_vb) {

  z_vb <- 1 / exp(-(logP + log_Phi_xi_vb - log_1P - log_1_Phi_xi_vb))
  return(z_vb)
}

update_rho_ <- function(X, xi_vb) {

  log_pnorm <- pnorm(xi_vb, log.p = TRUE)
  log_1_pnorm <- pnorm(xi_vb, log.p = TRUE, lower.tail = FALSE)

  imr0 <- inv_mills_ratio_(0, xi_vb, log_1_pnorm, log_pnorm)
  rho_vb <- z_vb * (inv_mills_ratio_(1, xi_vb, log_1_pnorm, log_pnorm) - imr0) + imr0 + xi_vb

  return(rho_vb)
}

update_sig2_gam_0_ <- function(sig_0_inv2, d) {

  sig2_gam_0_vb <- rep(1 / (1 + sig_0_inv2), d)

  return(sig2_gam_0_vb)
}

update_mu_gam_0_ <- function(sig_0_inv2, eta_0, X, mu_gam_1_vb, rho_vb) {

  mu_gam_0_vb <- 1/sig2_gam_0_vb * (eta_0 - colSums(sweep(X, 2, mu_gam_1_vb, "*") - rho_vb)) # X is n x d, mu_gam_1_vb is vector of length d, rho_vb is n x d

  return(mu_gam_0_vb)
}

update_sig2_gam_1_ <- function(X, tau_inv2_vb, lambda_inv2_vb) {

  sig2_gam_1_vb <- 1 / (colSums(X^2) + tau_inv2_vb * lambda_inv2_vb)

  return(sig2_gam_1_vb)
}

update_mu_gam_1_ <- function(sig2_gam_1_vb, X, rho_vb, mu_gam_0_vb) {

  mu_gam_1_vb <- 1/sig2_gam_1_vb * colSums(X * sweep(rho_vb, 2, mu_gam_0_vb, "-"))

  return(mu_gam_1_vb)
}

update_eta_lambda_ <- function(tau_inv2_vb, sig2_gam_1_vb, mu_gam_1_vb, b_inv_vb) {

  eta_lambda <- 1/2 * tau_inv2_vb * (sig2_gam_1_vb + mu_gam_1_vb) + b_inv_vb # these are just vectors, so no need for sweep() just yet

  return(eta_lambda)
}

update_lambda_inv2_ <- function(eta_lambda, log) {

  nu_lambda <- 1/2

  if(log) {

    lambda_inv2_vb <- nu_lambda / eta_lambda

    return(lambda_inv2_vb)
  }
  else {

    log_lambda_inv2_vb <- digamma(nu_lambda) - log(eta_lambda)

    return(log_lambda_inv2_vb)
  }
}

update_eta_tau_ <- function(lambda_inv2_vb, sig2_gam_1_vb, mu_gam_1_vb, a_inv_vb) {

  eta_tau <- 1/2 * sum(lambda_inv2_vb * (sig2_gam_1_vb + mu_gam_1_vb)) + a_inv_vb

  return(eta_tau)
}

update_tau_inv2_ <- function(eta_tau, log) {

  nu_tau <- 1/2

  if(log) {

    tau_inv2_vb <- nu_tau / eta_tau

    return(tau_inv2_vb)
  }
  else {

    log_tau_inv2_vb <- digamma(nu_tau) / log(eta_tau)

    return(log_tau_inv2_vb)
  }
}

update_a_inv_ <- function(tau_inv2_vb, log) {

  c_a <- 1
  d_a <- tau_inv2_vb + 1

  if(log) {

    a_inv_vb <- c_a / d_a

    return(a_inv_vb)
  }
  else {

    log_a_inv_vb <- digamma(c_a) - log(d_a)

    return(log_a_inv_vb)
  }
}

update_b_inv_ <- function(lambda_inv2_vb, d, log) {

  c_b <- rep(1, d)
  d_b <- lambda_inv2_vb + 1

  if(log) {

    b_inv_vb <- c_b / d_b

    return(b_inv_vb)
  }
  else {

    log_b_inv_vb <- digamma(c_b) - log(d_b)

    return(log_b_inv_vb)
  }
}


