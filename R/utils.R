# Beta densities but infinities are replaced with 0's
# dbeta_no_NAs_ <- function(X, alpha, beta) {
#
#   db <- dbeta(X, alpha, beta)
#   db <- ifelse(is.infinite(db), 0, db)
# }

# Enforce labels on \mu_0, \mu in fit_mixture_model_()
enforce_labels_ <- function(shape_params) {

  d <- nrow(shape_params)

  alpha_0 <- shape_params[, 1]
  beta_0 <- shape_params[, 2]
  alpha_1 <- shape_params[, 3]
  beta_1 <- shape_params[, 4]

  swap_params_ <- function(a_0, b_0, a_1, b_1) {

    m_t_0 <- a_0 / (a_0 + b_0)
    m_t_1 <- a_1 / (a_1 + b_1)

    if(m_t_0 > m_t_1) {

      temp_a <- a_0
      a_0 <- a_1
      a_1 <- temp_a

      temp_b <- b_0
      b_0 <- b_1
      b_1 <- temp_b

    }
    return(c(a_0, b_0, a_1, b_1))
  }

  labeled_shape_params <- mapply(FUN = swap_params_,
                                 a_0 = alpha_0,
                                 b_0 = beta_0,
                                 a_1 = alpha_1,
                                 b_1 = beta_1)

  labeled_shape_params <- unlist(labeled_shape_params)
  labeled_shape_params <- t(matrix(labeled_shape_params, nrow = 4, ncol = d))

  return(labeled_shape_params)
}
