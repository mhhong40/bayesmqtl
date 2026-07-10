# Implemented in the full multivariate regression setting (p, d \geq 1)

update_z_ <- function(logP_upper, logP_lower, log_Phi_xi_vb, log_1_Phi_xi_vb) {

  z_vb <- plogis(logP_upper + log_Phi_xi_vb - logP_lower - log_1_Phi_xi_vb) # plogis() for added numerical stability
  return(z_vb)
}

update_rho_ <- function(X, z_vb, xi_vb) {

  imr0 <- inv_mills_ratio_(0, xi_vb)
  rho_vb <- z_vb * (inv_mills_ratio_(1, xi_vb) - imr0) + imr0 + xi_vb

  return(rho_vb)
}

update_sig2_gam_0_ <- function(sig_0_inv2, n, d) {

  sig2_gam_0_vb <- rep(1 / (n + sig_0_inv2), d)

  return(sig2_gam_0_vb)
}

update_mu_gam_0_ <- function(sig_0_inv2, sig2_gam_0_vb, eta_0, X, mu_gam_1_vb, rho_vb) {

  mu_gam_0_vb <- sig2_gam_0_vb * (sig_0_inv2 * eta_0 + colSums(rho_vb - X %*% mu_gam_1_vb))

  return(mu_gam_0_vb)
}

update_sig2_gam_1_ <- function(X, tau_inv2_vb, lambda_inv2_vb) {

  sig2_gam_1_vb <- 1 / (sweep(tcrossprod(lambda_inv2_vb, tau_inv2_vb), 1, colSums(X^2), "+"))

  return(sig2_gam_1_vb)
}

update_mu_gam_1_ <- function(sig2_gam_1_vb, X, d, rho_vb, mu_gam_0_vb, mu_gam_1_vb, resid_mu_gam_1, vectorized = TRUE) {

  p <- ncol(X)

  rho_min_all_mu <- sweep(rho_vb, 2, mu_gam_0_vb, "-") - X %*% mu_gam_1_vb

  if(vectorized) {

    for (j in 1:p) {

      mu_gam_1_j <- mu_gam_1_vb[j, ]
      mu_gam_1_vb[j, ] <- sig2_gam_1_vb[j, ] * crossprod(X[, j], rho_min_all_mu + outer(X[, j], mu_gam_1_j))
      rho_min_all_mu <- rho_min_all_mu - outer(X[, j], mu_gam_1_vb[j, ] - mu_gam_1_j)

    }
  }
  else {

    for (k in 1:d) { # Slower

      for (j in 1:p) {

        resid_mu_gam_1[, k] <- rho_vb[, k] - mu_gam_0_vb[k] - X %*% mu_gam_1_vb[, k] + mu_gam_1_vb[j, k] * X[, j, drop = FALSE]

        mu_gam_1_vb[j, k] <- sig2_gam_1_vb[j, k] * sum(X[, j] * resid_mu_gam_1[, k])
      }
    }
  }

  return(mu_gam_1_vb)
}

update_eta_lambda_ <- function(tau_inv2_vb, sig2_gam_1_vb, mu_gam_1_vb, a_inv_vb) {

  eta_lambda <- 1/2 * rowSums(sweep(sig2_gam_1_vb + mu_gam_1_vb^2, 2, tau_inv2_vb, "*")) + a_inv_vb

  return(eta_lambda)
}

update_lambda_inv2_ <- function(nu_lambda, eta_lambda, log) {

  if(log) {

    log_lambda_inv2_vb <- digamma(nu_lambda) - log(eta_lambda)

    return(log_lambda_inv2_vb)
  }
  else {

    lambda_inv2_vb <- nu_lambda / eta_lambda

    return(lambda_inv2_vb)
  }
}

update_eta_tau_ <- function(lambda_inv2_vb, sig2_gam_1_vb, mu_gam_1_vb, b_inv_vb) {

  eta_tau <- 1/2 * colSums(sweep(sig2_gam_1_vb + mu_gam_1_vb^2, 1, lambda_inv2_vb, "*")) + b_inv_vb

  return(eta_tau)
}

update_tau_inv2_ <- function(nu_tau, eta_tau, log) {

  if(log) {

    log_tau_inv2_vb <- digamma(nu_tau) - log(eta_tau)

    return(log_tau_inv2_vb)
  }
  else {

    tau_inv2_vb <- nu_tau / eta_tau

    return(tau_inv2_vb)
  }
}

update_a_inv_ <- function(c_a, d_a, log) {

  if(log) {

    log_a_inv_vb <- digamma(c_a) - log(d_a)

    return(log_a_inv_vb)
  }
  else {

    a_inv_vb <- c_a / d_a

    return(a_inv_vb)
  }
}

update_b_inv_ <- function(c_b, d_b, log) {

  if(log) {

    log_b_inv_vb <- digamma(c_b) - log(d_b)

    return(log_b_inv_vb)
  }
  else {

    b_inv_vb <- c_b / d_b

    return(b_inv_vb)
  }
}


