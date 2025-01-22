# Internal functions for each parameter update written here
# for abstraction in higher-level wrappers/readability.

update_m2_theta_vb_ <- function(mu_theta_vb, sig2_theta_vb) {

  mu_theta_vb^2 + sig2_theta_vb
}

update_sig2_theta_vb_ <- function(n, tau_vb, sig2_vb, c = 1) { # should be matrix!

  1 / (c * sweep(sig2_vb, 2, (n - 1)*tau_vb, "+"))

}

## tau's updates
update_eta_vb_ <- function(n, eta, c = 1) c * (n/2 + eta - 1) + 1

update_kappa_vb_ <- function(z, kappa, mat_x_m1, m2_theta, theta_vb, c = 1) {

  n <- nrow(z)

  c * (kappa + ( colSums(z^2) - 2 * colSums(z * mat_x_m1)  +
                   colSums(mat_x_m1^2) + (n - 1) * colSums(m2_theta) -
                   (n - 1) * colSums(theta_vb^2) )/ 2)
}

update_log_tau_vb_ <- function(eta_vb, kappa_vb) digamma(eta_vb) - log(kappa_vb)

## lambda's update
update_lambda_vb_ <- function(lambda, m2_theta, c = 1) {

  c * (1/2 * rowSums(m2_theta) + lambda)
}


