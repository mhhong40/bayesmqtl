#' Internal functions calculating the remaining ELBO terms for elbo() in bayesmqtl_core script.

# E log p(z | rest)
e_z_ <- function(n, kappa, kappa_vb, log_tau_vb, tau_vb, m2_theta) {

  sum( -n / 2 * log(2 * pi) + n / 2 * log_tau_vb - tau_vb *
    (kappa_vb - 1/2*colSums(m2_theta) - kappa) )

}

# E log p(tau | rest) - E log q(tau)
e_tau_ <- function(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb) {

  eps <- .Machine$double.eps^0.5

  sum( (eta - eta_vb) * log_tau_vb - (kappa - kappa_vb) * tau_vb +
        eta * log(kappa + eps) - eta_vb * log(kappa_vb + eps) - lgamma(eta) + lgamma(eta_vb) )
}

# E log p(theta | rest) - E log q(theta)
e_theta_ <- function(sig2_theta_vb, log_sig2_theta_vb, m2_theta) {

  # note that log_sig2_theta_vb is approx. = log(sig2_theta_vb) - 1/2
  1/2 * sum( -log(2 * pi) - log_sig2_theta_vb + (sig2_theta_vb * m2_theta)/2 + 1 )

}

# E log p(sig2_inv | rest) - E log q(sig2_inv)
e_sig2_ <- function(lambda, lambda_vb, sig2_vb) {

  sum( lambda - lambda_vb + (lambda_vb - lambda)*sig2_vb )
}
