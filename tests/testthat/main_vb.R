library(maxLik)
devtools::load_all()

# CpG data
dat <- readRDS("easy_simulation_datasets.RDS")
dat <- dat[[1]]

# SNP data
library(echoseq)
snps <- generate_snps(n = 250, p = 10, user_seed = NULL)
snps <- snps$snps

list_hyper <- auto_set_hyper_()
list_init <- auto_set_init_(Y = dat, X = snps, user_seed = NULL)

devtools::load_all()
bayesmqtl <- bayesmqtl(Y = dat, X = snps, list_hyper, list_init, tol_mix = 0.1, tol_vb = 0.1, maxit = 1000, verbose = FALSE)
