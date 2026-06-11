# Initialize parameters according to user input
# set_init_ <- function() {}

# Automatically initialize parameters
auto_set_init_ <- function(Y, X, user_seed) {

  p <- ncol(X)
  d <- ncol(Y)

  set.seed(user_seed)

  mu_gam_0_vb <- rnorm(d, mean = 0, sd = 1)
  mu_gam_1_vb <- matrix(rnorm(p*d, mean = 0, sd = 1/sqrt(p)), nrow = p, ncol = d, dimnames = list(colnames(X), colnames(Y)))

  tau_inv2_vb <- rgamma(d, shape = max(p, d), rate = 1)
  lambda_inv2_vb <- rgamma(p, shape = max(p, d), rate = 1) # Might have to revise, but keeping it the same as tau_inv2_vb for now since both have the same inv. Gamma parametrization

  list_init <- create_named_list_(mu_gam_0_vb, mu_gam_1_vb, tau_inv2_vb, lambda_inv2_vb)

  return(list_init)
}

# Set hyperparameters according to user input
# set_hyper_ <- function() {}

# Automatically set hyperparameters
auto_set_hyper_ <- function() {

  eta_0 <- 0
  sig_0_inv2 <- 0.5

  list_hyper <- create_named_list_(eta_0, sig_0_inv2)

  return(list_hyper)
}
