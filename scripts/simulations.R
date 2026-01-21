### Simulated data generation/testing script

## A number of general functions for data generation ----
generate_phi_ <- function(b, mu, mult) { # b is a scalar; mu and mult can be vector inputs

  del <- 1e-4 # An arbitrary choice; this can be tweaked
  phi_min <- del # Ensures that the distribution is not degenerate

  ind_l <- which(mu < 0.5)
  ind_u <- which(mu >= 0.5)

  phi_max <- rep(0, 5)

  phi_max[ind_l] <- mu[ind_l]^2 * ((1 - mu[ind_l]) / (1 + mu[ind_l])) # Ensures that \beta > \alpha > 1 if \mu < 0.5
  phi_max[ind_u] <- mu[ind_u] * ((1 - mu[ind_u])^2 / (2 - mu[ind_u])) # Ensures that \alpha > \beta > 1 if \mu >= 0.5

  # The input variable mult is used to further restrict the final variance of the simulated dataset
  phi <- runif(b, min = phi_min, max = pmax(phi_max*mult, del)) # b = batch size
  return(phi)
}

compute_alpha_ <- function(mu, phi){

  return(mu * (mu * (1 - mu) / phi - 1))
}

compute_beta_ <- function(mu, phi){

  return((1 - mu) * (mu * (1 - mu) / phi - 1))
}

# This operates on individual CpGs, meaning n is should be a scalar and par a single vector
rbeta_mixture_ <- function(n, par) {

  pi <- par[1]
  alpha_0 <- par[2]
  beta_0 <- par[3]
  alpha_1 <- par[4]
  beta_1 <- par[5]

  n_lower <- round((1 - pi)*n)
  n_upper <- n - n_lower

  return(c(rbeta(n_lower, alpha_0, beta_0), rbeta(n_upper, alpha_1, beta_1)))
}

## A number of general functions for evaluating accuracy and/or runtime ----
# Merge results from each CpG into one matrix
merge_results_ <- function(results_list) {

  results_mat <- lapply(results_list, FUN = function(x) {

    final_it <- length(x)
    return(x[[final_it]]$param_estimates)
  })  # Returns a list

  # Convert to matrix
  results_mat <- do.call(rbind, results_mat)
  return(results_mat[, -ncol(results_mat)]) # Doesn't include the current_LL column in final output
}

calculate_moments_ <- function(par_mat) {

  alpha_0 <- par_mat[, 2]
  beta_0 <- par_mat[, 3]
  alpha_1 <- par_mat[, 4]
  beta_1 <- par_mat[, 5]

  mu_0 <- alpha_0 / (alpha_0 + beta_0)
  phi_0 <- (alpha_0 * beta_0) / ((alpha_0 * beta_0)^2 * (alpha_0 + beta_0 + 1))
  mu_1 <- alpha_1 / (alpha_1 + beta_1)
  phi_1 <- (alpha_1 * beta_1) / ((alpha_1 * beta_1)^2 * (alpha_1 + beta_1 + 1))

  mom_mat <- cbind(par_mat[, 1], mu_0, phi_0, mu_1, phi_0)
  colnames(mom_mat) <- c("pi", "mu_0", "phi_0", "mu_1", "phi_0")
  return(mom_mat)
}

# General function to obtain CpG classifications
# To obtain true classifications (the so-called "key"), use par_mat from data generation;
# for predicted, use results
post_probs_ <- function(par_mat, Y) {

  pi_t <- par_mat[, 1]
  alpha_0_t <- par_mat[, 2]
  beta_0_t <- par_mat[, 3]
  alpha_1_t <- par_mat[, 4]
  beta_1_t <- par_mat[, 5]

  dens_lower <- dbeta(Y, alpha_0_t, beta_0_t, log = FALSE)
  dens_upper <- dbeta(Y, alpha_1_t, beta_1_t, log = FALSE)

  # Compute posterior P(observation belongs to higher mixture component) for all combinations of observations/CpGs
  post <- sweep(dens_upper, 1, pi_t, '*') / (sweep(dens_lower, 1, (1 - pi_t), '*') + sweep(dens_upper, 1, pi_t, '*'))
  return(post)
}

# Collect time elapsed from each iteration of the algorithm for each CpG
extract_runtimes_ <- function(results_list) {

  # Runtime data for a given results matrix is a 5-component list
  # Each component in the list is a vector of iteration runtimes
  # Number of iterations until convergence is simply given by (length of the vector + 1)
  times_list <- lapply(results_list, FUN = function(x) {

    return(sapply(x, function(y) { as.numeric(y[[2]]) }))
  })  # Returns a list

  times_list <- lapply(times_list, FUN = function(x) { x <- x[-1] }) # First iteration doesn't have a corresponding runtime recorded

  return(times_list)
}

# Collect log-likelihoods for a single CpG
extract_LLs_ <- function(single_results_list) {

  LL_list <- lapply(single_results_list, FUN = function(x) {

    if(length(x) != 2) {  # Handles the fact that the first element comprises param_estimates only (and thus has a length of 6)

      return(x[6])
    }
    else {

      LL <- x$param_estimates
      return(LL[6]) # current_LL is given by the final element of param_estimates
    }
  })  # Returns a list

  return(unlist(LL_list)) # Return vector output
}

total_runtimes_ <- function(runtime_list) {

  return(unlist(lapply(runtime_list, sum)))
}

total_iterations_ <- function(runtime_list) {

  return(unlist(lapply(runtime_list, length)))
}

## Actual data generation begins here ----
set.seed(12345) # For reproducibility

# Set d = number of CpGs, b = batch size (for each difficulty level)
d = 15
b = d/3

# Choose \pi based on easy/medium/hard settings
pi_easy <- runif(b, min = 0.4, max = 0.5)
pi_med <- runif(b, min = 0.2, max = 0.4)
pi_hard <- runif(b, min = 0.01, max = 0.2)

# Choose means for each of d = 15 CpGs
eps <- 1e-4

mu_0_easy <- runif(b, min = eps, max = 0.5)
min_easy <- pmin(0.3, (1 - eps) - mu_0_easy) # 0.3 as a bound ensures that none of the means are too close to 1, but it can probably be slightly increased
max_easy <- pmax(0.3, (1 - eps) - mu_0_easy)
mu_1_easy <- mu_0_easy + runif(b, min = min_easy, max = max_easy) # Enforces \mu_0 < \mu_1

mu_0_med <- runif(b, min = eps, max = 0.75)
min_med <- pmin(0.25, 0.5)
mu_1_med <- mu_0_med + runif(b, min = 0.25, max = 0.5)

mu_0_hard <- runif(b, min = eps, max = 0.8)
mu_1_hard <- mu_0_hard + runif(b, min = 0.01, max = 0.20)

# Choose variances for each of d = 15 CpGs
phi_0_easy <- generate_phi_(b, mu_0_easy, mult = runif(b, min = 0.25, max = 0.5))
phi_1_easy <- generate_phi_(b, mu_1_easy, mult = runif(b, min = 0.25, max = 0.5))

phi_0_med <- generate_phi_(b, mu_0_med, mult = runif(b, min = 0.6, max = 0.85))
phi_1_med <- generate_phi_(b, mu_1_med, mult = runif(b, min = 0.6, max = 0.85))

phi_0_hard <- generate_phi_(b, mu_0_hard, 1) # No additional restriction on variance
phi_1_hard <- generate_phi_(b, mu_1_hard, 1)

# Save moments just in case
mom_easy <- cbind(mu_0_easy, phi_0_easy, mu_1_easy, phi_1_easy)
saveRDS(mom_easy, file = "easy_simulation_moments.RDS")

mom_med <- cbind(mu_0_med, phi_0_med, mu_1_med, phi_1_med)
saveRDS(mom_med, file = "med_simulation_moments.RDS")

mom_hard <- cbind(mu_0_hard, phi_0_hard, mu_1_hard, phi_1_hard)
saveRDS(mom_hard, file = "hard_simulation_moments.RDS")

# Compute alphas and betas for each
alpha_0_easy <- compute_alpha_(mu_0_easy, phi_0_easy)
beta_0_easy <- compute_beta_(mu_0_easy, phi_0_easy)
alpha_1_easy <- compute_alpha_(mu_1_easy, phi_1_easy)
beta_1_easy <- compute_beta_(mu_1_easy, phi_1_easy)
par_easy <- cbind(pi_easy, alpha_0_easy, beta_0_easy, alpha_1_easy, beta_1_easy)
saveRDS(par_easy, file = "easy_simulation_params.RDS")

alpha_0_med <- compute_alpha_(mu_0_med, phi_0_med)
beta_0_med <- compute_beta_(mu_0_med, phi_0_med)
alpha_1_med <- compute_alpha_(mu_1_med, phi_1_med)
beta_1_med <- compute_beta_(mu_1_med, phi_1_med)
par_med <- cbind(pi_med, alpha_0_med, beta_0_med, alpha_1_med, beta_1_med)
saveRDS(par_med, file = "med_simulation_params.RDS")

alpha_0_hard <- compute_alpha_(mu_0_hard, phi_0_hard)
beta_0_hard <- compute_beta_(mu_0_hard, phi_0_hard)
alpha_1_hard <- compute_alpha_(mu_1_hard, phi_1_hard)
beta_1_hard <- compute_beta_(mu_1_hard, phi_1_hard)
par_hard <- cbind(pi_hard, alpha_0_hard, beta_0_hard, alpha_1_hard, beta_1_hard)
saveRDS(par_hard, file = "hard_simulation_params.RDS")

# Easy datasets
# par_easy <- readRDS("easy_simulation_params.RDS")

data_250_easy <- mapply(rbeta_mixture_,
                        n = rep(250, 5),
                        par = split(par_easy, f = row(par_easy)))

data_500_easy <- mapply(rbeta_mixture_,
                        n = rep(500, 5),
                        par = split(par_easy, f = row(par_easy)))

data_750_easy <- mapply(rbeta_mixture_,
                        n = rep(750, 5),
                        par = split(par_easy, f = row(par_easy)))

data_1000_easy <- mapply(rbeta_mixture_,
                         n = rep(1000, 5),
                         par = split(par_easy, f = row(par_easy)))

saveRDS(list(data_250_easy, data_500_easy, data_750_easy, data_1000_easy), file = "easy_simulation_datasets.RDS")

# Medium datasets
# par_med <- readRDS("med_simulation_params.RDS")

data_250_med <- mapply(rbeta_mixture_,
                        n = rep(250, 5),
                        par = split(par_med, f = row(par_med)))

data_500_med <- mapply(rbeta_mixture_,
                        n = rep(500, 5),
                        par = split(par_med, f = row(par_med)))

data_750_med <- mapply(rbeta_mixture_,
                        n = rep(750, 5),
                        par = split(par_med, f = row(par_med)))

data_1000_med <- mapply(rbeta_mixture_,
                         n = rep(1000, 5),
                         par = split(par_med, f = row(par_med)))

saveRDS(list(data_250_med, data_500_med, data_750_med, data_1000_med), file = "med_simulation_datasets.RDS")

# Hard datasets
# par_hard <- readRDS("hard_simulation_params.RDS")

data_250_hard <- mapply(rbeta_mixture_,
                       n = rep(250, 5),
                       par = split(par_hard, f = row(par_hard)))

data_500_hard <- mapply(rbeta_mixture_,
                       n = rep(500, 5),
                       par = split(par_hard, f = row(par_hard)))

data_750_hard <- mapply(rbeta_mixture_,
                       n = rep(750, 5),
                       par = split(par_hard, f = row(par_hard)))

data_1000_hard <- mapply(rbeta_mixture_,
                        n = rep(1000, 5),
                        par = split(par_hard, f = row(par_hard)))

saveRDS(list(data_250_hard, data_500_hard, data_750_hard, data_1000_hard), file = "hard_simulation_datasets.RDS")


## Running fit_mixture_model_() for each dataset ----
library(maxLik)
library(here)
library(easybio) # To use split_matrix()
devtools::load_all() # Loads bayesmqtl

## Easy data ----
easy_params <- readRDS("easy_simulation_params.RDS")
easy_data <- readRDS("easy_simulation_datasets.RDS")
# n = 250
e_d250_results_0.1 <- mapply(fit_mixture_model_,
                  Y = split_matrix(easy_data[[1]], chunk_size = 1),
                  parametrization = rep("shape", 5),
                  tol = rep(0.1, 5))

e_d250_results_0.01 <- mapply(fit_mixture_model_,
                              Y = split_matrix(easy_data[[1]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.01, 5))

e_d250_results_0.001 <- mapply(fit_mixture_model_,
                              Y = split_matrix(easy_data[[1]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.001, 5))

e_d250_results_0.0001 <- mapply(fit_mixture_model_,
                               Y = split_matrix(easy_data[[1]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.0001, 5))

# n = 500
e_d500_results_0.1 <- mapply(fit_mixture_model_,
                             Y = split_matrix(easy_data[[2]], chunk_size = 1),
                             parametrization = rep("shape", 5),
                             tol = rep(0.1, 5))

e_d500_results_0.01 <- mapply(fit_mixture_model_,
                              Y = split_matrix(easy_data[[2]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.01, 5))

e_d500_results_0.001 <- mapply(fit_mixture_model_,
                               Y = split_matrix(easy_data[[2]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.001, 5))

e_d500_results_0.0001 <- mapply(fit_mixture_model_,
                                Y = split_matrix(easy_data[[2]], chunk_size = 1),
                                parametrization = rep("shape", 5),
                                tol = rep(0.0001, 5))

# n = 750
e_d750_results_0.1 <- mapply(fit_mixture_model_,
                             Y = split_matrix(easy_data[[3]], chunk_size = 1),
                             parametrization = rep("shape", 5),
                             tol = rep(0.1, 5))

e_d750_results_0.01 <- mapply(fit_mixture_model_,
                              Y = split_matrix(easy_data[[3]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.01, 5))

e_d750_results_0.001 <- mapply(fit_mixture_model_,
                               Y = split_matrix(easy_data[[3]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.001, 5))

e_d750_results_0.0001 <- mapply(fit_mixture_model_,
                                Y = split_matrix(easy_data[[3]], chunk_size = 1),
                                parametrization = rep("shape", 5),
                                tol = rep(0.0001, 5))

# n = 1000
e_d1000_results_0.1 <- mapply(fit_mixture_model_,
                             Y = split_matrix(easy_data[[4]], chunk_size = 1),
                             parametrization = rep("shape", 5),
                             tol = rep(0.1, 5))

e_d1000_results_0.01 <- mapply(fit_mixture_model_,
                              Y = split_matrix(easy_data[[4]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.01, 5))

e_d1000_results_0.001 <- mapply(fit_mixture_model_,
                               Y = split_matrix(easy_data[[4]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.001, 5))

e_d1000_results_0.0001 <- mapply(fit_mixture_model_,
                                Y = split_matrix(easy_data[[4]], chunk_size = 1),
                                parametrization = rep("shape", 5),
                                tol = rep(0.0001, 5))
save(list = ls(all.names = TRUE), file = "easy_simulation_results.RDA")
rm(list = ls(all.names = TRUE)) # Clear the environment before running each set of simulations corresponding to a certain difficulty level

## Medium data ----
med_data <- readRDS("med_simulation_datasets.RDS")
# n = 250
m_d250_results_0.1 <- mapply(fit_mixture_model_,
                             Y = split_matrix(med_data[[1]], chunk_size = 1),
                             parametrization = rep("shape", 5),
                             tol = rep(0.1, 5))

m_d250_results_0.01 <- mapply(fit_mixture_model_,
                              Y = split_matrix(med_data[[1]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.01, 5))

m_d250_results_0.001 <- mapply(fit_mixture_model_,
                               Y = split_matrix(med_data[[1]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.001, 5))

m_d250_results_0.0001 <- mapply(fit_mixture_model_,
                                Y = split_matrix(med_data[[1]], chunk_size = 1),
                                parametrization = rep("shape", 5),
                                tol = rep(0.0001, 5))

# n = 500
m_d500_results_0.1 <- mapply(fit_mixture_model_,
                             Y = split_matrix(med_data[[2]], chunk_size = 1),
                             parametrization = rep("shape", 5),
                             tol = rep(0.1, 5))

m_d500_results_0.01 <- mapply(fit_mixture_model_,
                              Y = split_matrix(med_data[[2]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.01, 5))

m_d500_results_0.001 <- mapply(fit_mixture_model_,
                               Y = split_matrix(med_data[[2]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.001, 5))

m_d500_results_0.0001 <- mapply(fit_mixture_model_,
                                Y = split_matrix(med_data[[2]], chunk_size = 1),
                                parametrization = rep("shape", 5),
                                tol = rep(0.0001, 5)) # Hit maxit = 1000

# n = 750
m_d750_results_0.1 <- mapply(fit_mixture_model_,
                             Y = split_matrix(med_data[[3]], chunk_size = 1),
                             parametrization = rep("shape", 5),
                             tol = rep(0.1, 5))

m_d750_results_0.01 <- mapply(fit_mixture_model_,
                              Y = split_matrix(med_data[[3]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.01, 5))

m_d750_results_0.001 <- mapply(fit_mixture_model_,
                               Y = split_matrix(med_data[[3]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.001, 5))

m_d750_results_0.0001 <- mapply(fit_mixture_model_,
                                Y = split_matrix(med_data[[3]], chunk_size = 1),
                                parametrization = rep("shape", 5),
                                tol = rep(0.0001, 5))

# n = 1000
m_d1000_results_0.1 <- mapply(fit_mixture_model_,
                              Y = split_matrix(med_data[[4]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.1, 5))

m_d1000_results_0.01 <- mapply(fit_mixture_model_,
                               Y = split_matrix(med_data[[4]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.01, 5))

m_d1000_results_0.001 <- mapply(fit_mixture_model_,
                                Y = split_matrix(med_data[[4]], chunk_size = 1),
                                parametrization = rep("shape", 5),
                                tol = rep(0.001, 5))

m_d1000_results_0.0001 <- mapply(fit_mixture_model_,
                                 Y = split_matrix(med_data[[4]], chunk_size = 1),
                                 parametrization = rep("shape", 5),
                                 tol = rep(0.0001, 5))

save(list = ls(all.names = TRUE), file = "med_simulation_results.RDA") # Just note that currently, m_d500_results_0.0001 is missing
rm(list = ls(all.names = TRUE)) # Clear the environment before running each set of simulations corresponding to a certain difficulty level

## Hard data ----
hard_data <- readRDS("hard_simulation_datasets.RDS")
# n = 250
h_d250_results_0.1 <- mapply(fit_mixture_model_,
                             Y = split_matrix(hard_data[[1]], chunk_size = 1),
                             parametrization = rep("shape", 5),
                             tol = rep(0.1, 5))

h_d250_results_0.01 <- mapply(fit_mixture_model_,
                              Y = split_matrix(hard_data[[1]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.01, 5))

h_d250_results_0.001 <- mapply(fit_mixture_model_,
                               Y = split_matrix(hard_data[[1]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.001, 5))

h_d250_results_0.0001 <- mapply(fit_mixture_model_,
                                Y = split_matrix(hard_data[[1]], chunk_size = 1),
                                parametrization = rep("shape", 5),
                                tol = rep(0.0001, 5))

# n = 500
h_d500_results_0.1 <- mapply(fit_mixture_model_,
                             Y = split_matrix(hard_data[[2]], chunk_size = 1),
                             parametrization = rep("shape", 5),
                             tol = rep(0.1, 5))

h_d500_results_0.01 <- mapply(fit_mixture_model_,
                              Y = split_matrix(hard_data[[2]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.01, 5))

h_d500_results_0.001 <- mapply(fit_mixture_model_,
                               Y = split_matrix(hard_data[[2]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.001, 5))

h_d500_results_0.0001 <- mapply(fit_mixture_model_,
                                Y = split_matrix(hard_data[[2]], chunk_size = 1),
                                parametrization = rep("shape", 5),
                                tol = rep(0.0001, 5))

# n = 750
h_d750_results_0.1 <- mapply(fit_mixture_model_,
                             Y = split_matrix(hard_data[[3]], chunk_size = 1),
                             parametrization = rep("shape", 5),
                             tol = rep(0.1, 5))

h_d750_results_0.01 <- mapply(fit_mixture_model_,
                              Y = split_matrix(hard_data[[3]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.01, 5))

h_d750_results_0.001 <- mapply(fit_mixture_model_,
                               Y = split_matrix(hard_data[[3]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.001, 5))

h_d750_results_0.0001 <- lapply(split_matrix(hard_data[[3]], chunk_size = 1),
                                FUN = function(x) {
                                  return(fit_mixture_model_(x, parametrization = "shape", tol = 0.0001)) })

# n = 1000
h_d1000_results_0.1 <- mapply(fit_mixture_model_,
                              Y = split_matrix(hard_data[[4]], chunk_size = 1),
                              parametrization = rep("shape", 5),
                              tol = rep(0.1, 5))

h_d1000_results_0.01 <- mapply(fit_mixture_model_,
                               Y = split_matrix(hard_data[[4]], chunk_size = 1),
                               parametrization = rep("shape", 5),
                               tol = rep(0.01, 5))

h_d1000_results_0.001 <- mapply(fit_mixture_model_,
                                Y = split_matrix(hard_data[[4]], chunk_size = 1),
                                parametrization = rep("shape", 5),
                                tol = rep(0.001, 5))

h_d1000_results_0.0001 <- mapply(fit_mixture_model_,
                                 Y = split_matrix(hard_data[[4]], chunk_size = 1),
                                 parametrization = rep("shape", 5),
                                 tol = rep(0.0001, 5))

save(list = ls(all.names = TRUE), file = "hard_simulation_results.RDA")
rm(list = ls(all.names = TRUE)) # Clear the environment before running each set of simulations corresponding to a certain difficulty level



## Evaluating accuracy ----
easy_params <- readRDS("easy_simulation_params.RDS")
med_params <- readRDS("med_simulation_params.RDS")
hard_params <- readRDS("hard_simulation_params.RDS")
# 1. Error between predicted versus true distribution moments (mean and variance)
# all_easy_results <- mget(ls(pattern = "^e_d")) ## Easy
all_easy_p_est <- lapply(all_easy_results, merge_results_) # Parameter estimates for each dataset with d = 5
all_easy_p_est <- do.call(rbind, all_easy_p_est) # All easy results' parameter estimates merged

easy_mom <- readRDS("easy_simulation_moments.RDS")
easy_mom <- cbind(easy_params[, 1], easy_mom)
easy_mom_blocks <- rep(1, 16) %x% easy_mom

easy_mom_est <- calculate_moments_(all_easy_p_est)
easy_err <- easy_mom_est - easy_mom_blocks

# all_med_results <- mget(ls(pattern = "^m_d")) ## Medium
all_med_p_est <- lapply(all_med_results, merge_results_)
all_med_p_est <- do.call(rbind, all_med_p_est)

med_mom <- readRDS("med_simulation_moments.RDS")
med_mom <- cbind(med_params[, 1], med_mom)
med_mom_blocks <- rep(1, 16) %x% med_mom

med_mom_est <- calculate_moments_(all_med_p_est)
med_err <- med_mom_est - med_mom_blocks

# all_hard_results <- mget(ls(pattern = "^h_d")) ## Hard
all_hard_p_est <- lapply(all_hard_results, merge_results_)
all_hard_p_est <- do.call(rbind, all_hard_p_est)

hard_mom <- readRDS("hard_simulation_moments.RDS")
hard_mom <- cbind(hard_params[, 1], hard_mom)
hard_mom_blocks <- rep(1, 16) %x% hard_mom

hard_mom_est <- calculate_moments_(all_hard_p_est)
hard_err <- hard_mom_est - hard_mom_blocks # Interested in seeing whether means/variances are typically under or overestimated, so no absolute value

# Figure 6 (error distributions across all easy/medium/hard datasets) ----
pdf("Figure6.pdf", width = 15, height = 9)
par(mfrow = c(3, 5))
plot(density(easy_err[, 1]), main = "All easy",
     xlab = "Error in pi")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(easy_err[, 2]), main = "All easy",
     ylab = NA,
     xlab = "Error in mu_0")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(easy_err[, 3]), main = "All easy",
     ylab = NA,
     xlab = "Error in phi_0")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(easy_err[, 4]), main = "All easy",
     ylab = NA,
     xlab = "Error in mu_1")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(easy_err[, 5]), main = "All easy",
     ylab = NA,
     xlab = "Error in phi_1")
abline(v = 0, col = "red", lty = 2, lwd = 1)

plot(density(med_err[, 1]), main = "All medium",
     xlab = "Error in pi")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(med_err[, 2]), main = "All medium",
     ylab = NA,
     xlab = "Error in mu_0")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(med_err[, 3]), main = "All medium",
     ylab = NA,
     xlab = "Error in phi_0")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(med_err[, 4]), main = "All medium",
     ylab = NA,
     xlab = "Error in mu_1")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(med_err[, 5]), main = "All medium",
     ylab = NA,
     xlab = "Error in phi_1")
abline(v = 0, col = "red", lty = 2, lwd = 1)

plot(density(hard_err[, 1]), main = "All hard",
     xlab = "Error in pi")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(hard_err[, 2]), main = "All hard",
     ylab = NA,
     xlab = "Error in mu_0")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(hard_err[, 3]), main = "All hard",
     ylab = NA,
     xlab = "Error in phi_0")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(hard_err[, 4]), main = "All hard",
     ylab = NA,
     xlab = "Error in mu_1")
abline(v = 0, col = "red", lty = 2, lwd = 1)
plot(density(hard_err[, 5]), main = "All hard",
     ylab = NA,
     xlab = "Error in phi_1")
abline(v = 0, col = "red", lty = 2, lwd = 1)
dev.off()


# Figure 7 (changes in error for one CpG with varying n) ----
# Use hard V2, tol = 0.0001
x_n <- c(1000, 250, 500, 750)
y_n <- hard_err[c(2, 22, 42, 62), ] # <- varying n

pdf("Figure7.pdf", width = 15, height = 3)
par(mfrow = c(1, 5))
plot(x_n, y_n[, 1], main = "Error in pi for hard V2\ngiven various n, tol = 0.0001",
     xlab = "n", ylab = "Error",
     pch = 16, cex = 1, col = "red")
plot(x_n, y_n[, 2], main = "Error in mu_0 for hard V2\ngiven various n, tol = 0.0001",
     xlab = "n", ylab = "Error",
     pch = 16, cex = 1, col = "red")
plot(x_n, y_n[, 3], main = "Error in phi_0 for hard V2\ngiven various n, tol = 0.0001",
     xlab = "n", ylab = "Error",
     pch = 16, cex = 1, col = "red")
plot(x_n, y_n[, 4], main = "Error in mu_1 for hard V2\ngiven various n, tol = 0.0001",
     xlab = "n", ylab = "Error",
     pch = 16, cex = 1, col = "red")
plot(x_n, y_n[, 5], main = "Error in phi_1 for hard V2\ngiven various n, tol = 0.0001",
     xlab = "n", ylab = "Error",
     pch = 16, cex = 1, col = "red")
dev.off()

# Figure 8 (error distributions for one CpG with varying tol) ----
x_tol <- c(0.0001, 0.001, 0.01, 0.1)
x_tol_log <- log(x_tol, base = 10)
y_err_tol <- hard_err[c(2, 7, 12, 17), ] # <- varying tol

pdf("Figure8.pdf", width = 15, height = 3)
par(mfrow = c(1, 5))
plot(x_tol_log, y_err_tol[, 1], main = "Error in pi for hard V2\nusing various tol, n = 1000",
     xlab = "log_10(tol)", ylab = "Error",
     pch = 16, cex = 1, col = "red")
plot(x_tol_log, y_err_tol[, 2], main = "Error in mu_0 for hard V2\nusing various tol, n = 1000",
     xlab = "log_10(tol)", ylab = "Error",
     pch = 16, cex = 1, col = "red")
plot(x_tol_log, y_err_tol[, 3], main = "Error in phi_0 for hard V2\nusing various tol, n = 1000",
     xlab = "log_10(tol)", ylab = "Error",
     pch = 16, cex = 1, col = "red")
plot(x_tol_log, y_err_tol[, 4], main = "Error in mu_1 for hard V2\nusing various tol, n = 1000",
     xlab = "log_10(tol)", ylab = "Error",
     pch = 16, cex = 1, col = "red")
plot(x_tol_log, y_err_tol[, 5], main = "Error in phi_1 for hard V2\nusing various tol, n = 1000",
     xlab = "log_10(tol)", ylab = "Error",
     pch = 16, cex = 1, col = "red")
dev.off()

# 2. ROC curves (based on # CpGs correctly classified as belonging to the higher or lower component) ----
library(pROC)
# Figure 9 (ROC curves for hard V2, n = 1000, based on various tol)
pi_V2 <- hard_params[2, 1]; n <- 1000
n_lower <- round((1 - pi_V2)*n)
n_upper <- n - n_lower
true_classified <- c(rep(0, n_lower), rep(1, n_upper)) # Because of how rbeta_mixture_() is specified, the first n_lower observations correspond to the lower component, and the remaining ones correspond to the higher component

Y_mat <- hard_data[[4]] # Data matrix for hard data, n = 1000

results_mat_0.1 <- merge_results_(h_d1000_results_0.1)
post_prob_V2_0.1 <- post_probs_(results_mat_0.1, Y_mat)
post_prob_V2_0.1 <- post_prob_V2_0.1[, 2]

results_mat_0.01 <- merge_results_(h_d1000_results_0.01)
post_prob_V2_0.01 <- post_probs_(results_mat_0.01, Y_mat)
post_prob_V2_0.01 <- post_prob_V2_0.01[, 2]

results_mat_0.001 <- merge_results_(h_d1000_results_0.001)
post_prob_V2_0.001 <- post_probs_(results_mat_0.001, Y_mat)
post_prob_V2_0.001 <- post_prob_V2_0.001[, 2]

results_mat_0.0001 <- merge_results_(h_d1000_results_0.0001)
post_prob_V2_0.0001 <- post_probs_(results_mat_0.0001, Y_mat)
post_prob_V2_0.0001 <- post_prob_V2_0.0001[, 2]

roc_0.1 <- plot.roc(true_classified, post_prob_V2_0.1, percent = TRUE,
                    direction = "<",
                    smooth = TRUE, smooth.method = "binormal")
roc_0.01 <- plot.roc(true_classified, post_prob_V2_0.01,
                     direction = "<",
                     smooth = TRUE, smooth.method = "binormal")
roc_0.001 <- plot.roc(true_classified, post_prob_V2_0.001,
                      direction = "<",
                      smooth = TRUE, smooth.method = "binormal")
roc_0.0001 <- plot.roc(true_classified, post_prob_V2_0.0001,
                       direction = "<",
                       smooth = TRUE, smooth.method = "binormal")

pdf("Figure9.pdf", width = 5, height = 5)
plot(roc_0.1, col = "green4",
     main = "ROC curves for hard V2\nusing various tol, n = 1000")
lines(roc_0.01, col = "magenta3")
lines(roc_0.001, col = "orange")
lines(roc_0.0001, col = "blue")
legend(x = 30, y = 20,
       legend = c("tol = 0.0001", "tol = 0.001", "tol = 0.01", "tol = 0.1"),
       col = c("blue", "orange", "magenta3", "green4"),
       lty = rep(1, 4), bty = "n", cex = 0.7)
dev.off()

# Figure 10 (ROC curves for hard V2, tol = 0.0001, given various n)
# Need true classifications for other sample sizes. Recycling data objects here out of laziness
n <- 250
n_lower <- round((1 - pi_V2)*n)
n_upper <- n - n_lower
true_classified_250 <- c(rep(0, n_lower), rep(1, n_upper))
Y_mat <- hard_data[[1]] # Data matrix for hard data, n = 250

results_mat_250 <- merge_results_(h_d250_results_0.0001)
post_prob_V2_250 <- post_probs_(results_mat_250, Y_mat)
post_prob_V2_250 <- post_prob_V2_250[, 2]

n <- 500
n_lower <- round((1 - pi_V2)*n)
n_upper <- n - n_lower
true_classified_500 <- c(rep(0, n_lower), rep(1, n_upper))
Y_mat <- hard_data[[2]] # Data matrix for hard data, n = 250

results_mat_500 <- merge_results_(h_d500_results_0.0001)
post_prob_V2_500 <- post_probs_(results_mat_500, Y_mat)
post_prob_V2_500 <- post_prob_V2_500[, 2]

n <- 750
n_lower <- round((1 - pi_V2)*n)
n_upper <- n - n_lower
true_classified_750 <- c(rep(0, n_lower), rep(1, n_upper))
Y_mat <- hard_data[[3]] # Data matrix for hard data, n = 250

results_mat_750 <- merge_results_(h_d750_results_0.0001)
post_prob_V2_750 <- post_probs_(results_mat_750, Y_mat)
post_prob_V2_750 <- post_prob_V2_750[, 2]

roc_250 <- plot.roc(true_classified_250, post_prob_V2_250, percent = TRUE,
                    direction = "<",
                    smooth = TRUE, smooth.method = "binormal")

roc_500 <- plot.roc(true_classified_500, post_prob_V2_500, percent = TRUE,
                    direction = "<",
                    smooth = TRUE, smooth.method = "binormal")

roc_750 <- plot.roc(true_classified_750, post_prob_V2_750, percent = TRUE,
                    direction = "<",
                    smooth = TRUE, smooth.method = "binormal")

pdf("Figure10.pdf", width = 5, height = 5)
plot(roc_250, col = "blue",
     main = "ROC curves for hard V2\ngiven various n, tol = 0.0001")
lines(roc_500, col = "orange")
lines(roc_750, col = "magenta3")
lines(roc_0.0001, col = "green4") # Same as what roc_1000 would be
legend(x = 30, y = 20,
       legend = c("n = 250", "n = 500", "n = 750", "n = 1000"),
       col = c("blue", "orange", "magenta3", "green4"),
       lty = rep(1, 4), bty = "n", cex = 0.7)
dev.off()

## Evaluating runtime ----
# Fig 1: Examine distribution of total runtimes across all datasets of varying difficulty (easy, medium, hard)
load("easy_simulation_results.RDA")
all_easy_results <- mget(ls(pattern = "^e_d"))
all_easy_runtimes <- lapply(all_easy_results, extract_runtimes_) # Runtimes in list form
all_easy_runtimes <- unlist(lapply(all_easy_runtimes, total_runtimes_)) # Obtain totals in vector form

load("med_simulation_results.RDA")
all_med_results <- mget(ls(pattern = "^m_d"))
all_med_runtimes <- lapply(all_med_results, extract_runtimes_) # Runtimes in list form
all_med_runtimes <- unlist(lapply(all_med_runtimes, total_runtimes_)) # Obtain totals in vector form

load("hard_simulation_results.RDA")
all_hard_results <- mget(ls(pattern = "^h_d"))
all_hard_runtimes <- lapply(all_hard_results, extract_runtimes_) # Runtimes in list form
all_hard_runtimes <- unlist(lapply(all_hard_runtimes, total_runtimes_)) # Obtain totals in vector form

pdf("Figure1.pdf", width = 10, height = 3)
par(mfrow = c(1, 3))
hist(all_easy_runtimes, breaks = 20,
     main = "All easy",
     xlab = NA)
hist(all_med_runtimes, breaks = 20,
     main = "All medium",
     xlab = "Total runtime (s)")
hist(all_hard_runtimes, breaks = 20,
     main = "All hard",
     xlab = NA)
dev.off()

# Fig 1.5: Examine distribution of total runtimes across all datasets of varying values of n
runtimes_250 <- c(unlist(all_easy_runtimes[21:40]), unlist(all_med_runtimes[21:40]), unlist(all_hard_runtimes[21:40]))
runtimes_500 <- c(unlist(all_easy_runtimes[41:60]), unlist(all_med_runtimes[41:60]), unlist(all_hard_runtimes[41:60]))
runtimes_750 <- c(unlist(all_easy_runtimes[61:80]), unlist(all_med_runtimes[61:80]), unlist(all_hard_runtimes[61:80]))
runtimes_1000 <- c(unlist(all_easy_runtimes[1:20]), unlist(all_med_runtimes[1:20]), unlist(all_hard_runtimes[1:20]))

pdf("Figure1.5.pdf", width = 7, height = 6)
par(mfrow = c(2, 2))
hist(runtimes_250, breaks = 20,
     main = "All n = 250",
     xlab = "Total runtime (s)")
hist(runtimes_500, breaks = 20,
     main = "All n = 500",
     xlab = "Total runtime (s)")
hist(runtimes_750, breaks = 20,
     main = "All n = 750",
     xlab = "Total runtime (s)")
hist(runtimes_1000, breaks = 20,
     main = "All n = 1000",
     xlab = "Total runtime (s)")
dev.off()

# Fig 2: Examine distribution of number of iterations until convergence across all datasets of varying difficulty (easy, medium, hard)
all_easy_it_nums <- unlist(lapply(all_easy_results, total_iterations_))
all_med_it_nums <- unlist(lapply(all_med_results, total_iterations_))
all_hard_it_nums <- unlist(lapply(all_hard_results, total_iterations_))

pdf("Figure2.pdf", width = 10, height = 3)
par(mfrow = c(1, 3))
hist(all_easy_it_nums, breaks = 20,
     main = "All easy",
     xlab = NA)
hist(all_med_it_nums, breaks = 20,
     main = "All medium",
     xlab = "Final number of iterations")
hist(all_hard_it_nums, breaks = 20,
     main = "All hard",
     xlab = NA)
dev.off()

# Fig 2.5: Do the same but with varied n
it_nums_250 <- c(all_easy_it_nums[21:40], all_med_it_nums[21:40], all_hard_it_nums[21:40])
it_nums_500 <- c(all_easy_it_nums[41:60], all_med_it_nums[41:60], all_hard_it_nums[41:60])
it_nums_750 <- c(all_easy_it_nums[61:80], all_med_it_nums[61:80], all_hard_it_nums[61:80])
it_nums_1000 <- c(all_easy_it_nums[1:20], all_med_it_nums[1:20], all_hard_it_nums[1:20])

pdf("Figure2.5.pdf", width = 7, height = 6)
par(mfrow = c(2, 2))
hist(it_nums_250, breaks = 20,
     main = "All n = 250",
     xlab = "Final number of iterations")
hist(it_nums_500, breaks = 20,
     main = "All n = 500",
     xlab = "Final number of iterations")
hist(it_nums_750, breaks = 20,
     main = "All n = 750",
     xlab = "Final number of iterations")
hist(it_nums_1000, breaks = 20,
     main = "All n = 1000",
     xlab = "Final number of iterations")
dev.off()

# Figure 3: Cumulative sum plots of runtime versus iteration number
# 3A: Converged within 10 to 99 iterations (63 in this case)
time_short <- extract_runtimes_(m_d500_results_0.01)
time_short <- time_short[[3]]

# 3B: Converged within 100 to 500 iterations (247 in this case)
time_sm <- extract_runtimes_(m_d250_results_0.0001)
time_sm <- time_sm[[2]]

# 3C: Converged within 501 to 800 iterations (684 in this case)
load("hard_simulation_results.RDA")
time_ml <- extract_runtimes_(h_d1000_results_0.0001)
time_ml <- time_ml[[4]]

# 3D: Longest to converge (maxit = 1000 iterations)
time_max <- extract_runtimes_(m_d500_results_0.0001)
time_max <- time_max[[2]]

# pdf(file = "my_plots.pdf", width = 8, height = 11) <- replace with file name
pdf(file = "Figure3.pdf", width = 10, height = 10)
par(mfrow = c(2, 2))
plot(cumsum(time_short), type = "p",
     main = "Converged within 63 iterations \n(medium V3, n = 500, tol = 0.01)",
     xlab = "Iteration number", ylab = "Runtime until convergence (s)", cex = 0.7)
plot(cumsum(time_sm), type = "p",
     main = "Converged within 247 iterations \n(medium V2, n = 250, tol = 0.0001)",
     xlab = "Iteration number", ylab = "Runtime until convergence (s)", cex = 0.7)
plot(cumsum(time_ml), type = "p",
     main = "Converged within 684 iterations \n(hard V4, n = 1000, tol = 0.0001)",
     xlab = "Iteration number", ylab = "Runtime until convergence (s)", cex = 0.7)
plot(cumsum(time_max), type = "p",
     main = "Not converged after maxit = 1000 iterations \n(medium V2, n = 500, tol = 0.0001)",
     xlab = "Iteration number", ylab = "Runtime until convergence (s)", cex = 0.7)
dev.off()

# Figure 4: Compare course of log-likelihood of one CpG between various values of tol
ll_0.1 <- extract_LLs_(h_d1000_results_0.1[[2]]) # Use V2 of hard dataset, n = 1000
ll_0.01 <- extract_LLs_(h_d1000_results_0.01[[2]])
ll_0.001 <- extract_LLs_(h_d1000_results_0.001[[2]])
ll_0.0001 <- extract_LLs_(h_d1000_results_0.0001[[2]])

x_4 <- c(length(ll_0.1), length(ll_0.01), length(ll_0.001), length(ll_0.0001))
y_4 <- c(last(ll_0.1), last(ll_0.01), last(ll_0.001), last(ll_0.0001))

pdf(file = "Figure4.pdf", width = 5, height = 5)
plot(ll_0.0001, type = "l",
     main = "Course of log-likelihood for hard V2\nusing various tol, n = 1000",
     xlab = "Iteration number",
     ylab = "Log-likelihood", ylim = c(675, 681),
     col = "blue")
lines(ll_0.001, col = "orange")
lines(ll_0.01, col = "magenta3")
lines(ll_0.1, col = "green4")
points(x_4, y_4, pch = 4, cex = 0.5)
legend(x = 450, y = 676.3,
       legend = c("tol = 0.0001", "tol = 0.001", "tol = 0.01", "tol = 0.1"),
       col = c("blue", "orange", "magenta3", "green4"),
       lty = rep(1, 4), bty = "n", cex = 0.7)
dev.off()

# Figure 5: Compare course of log-likelihood of one CpG between various values of n, use tol = 0.0001
library(dplyr)
ll_250 <- extract_LLs_(h_d250_results_0.0001[[2]])
ll_500 <- extract_LLs_(h_d500_results_0.0001[[2]])
ll_750 <- extract_LLs_(h_d750_results_0.0001[[2]])
ll_1000 <- extract_LLs_(h_d1000_results_0.0001[[2]])

x_5 <- c(length(ll_250), length(ll_500), length(ll_750), length(ll_1000))
y_5 <- c(last(ll_250), last(ll_500), last(ll_750), last(ll_1000))

pdf(file = "Figure5.pdf", width = 5, height = 5)
plot(ll_1000, type = "l",
     main = "Course of log-likelihood for hard V2\ngiven various n, tol = 0.0001",
     xlab = "Iteration number",
     ylab = "Log-likelihood", ylim = c(0, 681),
     col = "green4")
lines(ll_750, col = "magenta3")
lines(ll_500, col = "orange")
lines(ll_250, col = "blue")
points(x_5, y_5, pch = 4, cex = 0.5)
legend(x = 450, y = 150, fill = NULL, border = "white",
       legend = c("n = 250", "n = 500", "n = 750", "n = 1000"),
       col = c("blue", "orange", "magenta3", "green4"),
       lty = rep(1, 4), bty = "n", cex = 0.7)
dev.off()
