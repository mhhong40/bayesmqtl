library(maxLik)
devtools::load_all()

# CpG data for unassociated CpGs
dat <- readRDS("easy_simulation_datasets.RDS")
dat <- dat[[1]]
Y_unassoc <- dat[, c(1:3)]

# SNP data
library(echoseq)
snps <- generate_snps(n = 250, p = 15, user_seed = 5)
snps <- snps$snps

# Generate CpG data for associated CpGs
n <- nrow(snps)
p <- ncol(snps)
d <- 5
d_assoc <- 2

library(extraDistr) # For rbern()
set.seed(5)
true_mu_gam_1 <- matrix(rnorm(p*d_assoc), nrow = p, ncol = d_assoc) * matrix(rbern(p*d_assoc, prob = 0.3), nrow = p, ncol = d_assoc)
true_mu_gam_0 <- rnorm(d_assoc)

true_xi <- sweep(snps %*% true_mu_gam_1, 1, true_mu_gam_0, "+")

true_rho <- matrix(rnorm(n = n*d_assoc, mean = true_xi), nrow = n, ncol = d_assoc)
true_z <- ifelse(true_rho > 0, 1, 0) # True mixture component labels

shape_params <- readRDS("easy_simulation_params.RDS")
shape_params <- as.data.frame(shape_params[4:5, 2:5]) # Remove unassociated CpGs

Y_assoc_lower <- mapply(FUN = rbeta, n = rep(n, d_assoc), shape1 = shape_params$alpha_0_easy, shape2 = shape_params$beta_0_easy) # Lower component observations
Y_assoc_higher <- mapply(FUN = rbeta, n = rep(n, d_assoc), shape1 = shape_params$alpha_1_easy, shape2 = shape_params$beta_1_easy) # Higher

Y_assoc <- (1 - true_z)*Y_assoc_lower + true_z*Y_assoc_higher

Y <- cbind(Y_assoc, Y_unassoc)

# Run bayesmqtl
# Try just associated CpGs for now...
list_hyper <- auto_set_hyper_()
list_init <- auto_set_init_(Y = Y_assoc, X = snps, user_seed = 5)

devtools::load_all()
bayesmqtl <- bayesmqtl(Y = Y_assoc, X = snps, list_hyper, list_init, tol_mix = 0.1, tol_vb = 0.1, maxit = 1000, verbose = FALSE)
