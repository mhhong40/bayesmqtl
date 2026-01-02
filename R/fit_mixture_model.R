#' Contains the various functions needed to fit two-component Beta mixture models
#' to DNA methylation (microarray) data.
#'
#' Fits the mixture model via an expectation-maximization (EM) algorithm.
#' E and M-steps are both coded from scratch to facilitate parallel execution.
#' @param Y The methylation data matrix of size n x d, where n = the number of samples and d = the number of methylation sites.
#' @param threshold A value ranging from 0 to 1 (exclusive) indicating how the data distribution should first be split for parameter initialization. Default is 0.5 (data is split into upper and lower halves).
#' @param parameterization A string that is either "shape" or "md." Specifying "shape" performs the model fitting based on the alpha/beta shape parameterization of the Beta distribution, while md performs that based on the mean/dispersion parameterization.
#' @param ... More stuff here.
#' @param tol Tolerance value for the difference in the previous and current iterations' log-likelihood values that determines when the EM algorithm should stop. Default is 0.1.
#' @param maxit Maximum allowed number of iterations for the EM algorithm before it force quits.
fit_mixture_model_ <- function(Y, threshold = 0.5, parameterization = c("shape", "md"), tol = 0.1, maxit = 1000) {

  n <- nrow(Y)
  d <- ncol(Y)

  # Control over numerical under/overflow
  eps <- 1e-6
  Y <- pmin(pmax(Y, eps), 1 - eps)

  cutoffs <- matrixStats::colQuantiles(Y, probs = threshold)

  lower_mask <- sweep(Y, 2, cutoffs, '>=')
  lower <- Y * lower_mask
  upper <- Y * (1 - lower_mask)

  # Initialize model parameters (method of moments initialization)
  m_l <- colMeans(lower)
  v_l <- matrixStats::colVars(lower)

  m_u <- colMeans(upper)
  v_u <- matrixStats::colVars(upper)

  if (parameterization == "shape") {

    # Initialize model parameters (method of moments initialization for shape parameters)
    mix_params <- rep(threshold, d)

    alpha_0 <- m_l * (m_l * (1 - m_l) / v_l - 1)
    beta_0 <- (1 - m_l) * (m_l * (1 - m_l) / v_l - 1)

    alpha_1 <- m_u * (m_u * (1 - m_u) / v_u - 1)
    beta_1 <- (1 - m_u) * (m_u * (1 - m_u) / v_u - 1)

    current_LL <- colSums(bm_density_(Y, alpha_0, beta_0, alpha_1, beta_1, mix_params, log = TRUE))
    prev_LL <- current_LL - 100 # arbitrary starting point

    result <- cbind(mix_params, alpha_0, beta_0, alpha_1, beta_1, current_LL)

    iteration <- 1
    iteration_results <- list()
    iteration_results[[iteration]] <- result

    # Only run EM algorithm on CpGs that have not yet converged
    not_converged <- rep(TRUE, d)

    while(length(not_converged) > 0 | iteration > maxit) { # Runs the algorithm until all models have reached convergence

      iteration <- iteration + 1

      # E-step
      mix_params_current <- e_step_(Y[, not_converged, drop = FALSE], mix_params[not_converged], alpha_0[not_converged], beta_0[not_converged], alpha_1[not_converged], beta_1[not_converged])

      # M-step
      updated_shape_params <- m_step_(Y[, not_converged, drop = FALSE], mix_params[not_converged], alpha_0[not_converged], beta_0[not_converged], alpha_1[not_converged], beta_1[not_converged]) # need to create shape_params vector, accessing methods

      mix_params[not_converged] <- mix_params_current
      alpha_0[not_converged] <- updated_shape_params[, 1]
      beta_0[not_converged] <- updated_shape_params[, 2]
      alpha_1[not_converged] <- updated_shape_params[, 3]
      beta_1[not_converged] <- updated_shape_params[, 4]

      prev_LL <- current_LL
      current_LL <- colSums(bm_density_(Y, alpha_0, beta_0, alpha_1, beta_1, mix_params, log = TRUE)) # With updated parameters

      # Update the indices of CpGs that still have not converged
      not_converged <- abs(current_LL - prev_LL) > tol
      not_converged <- which(not_converged)

      # Storing step results
      result <- cbind(mix_params, alpha_0, beta_0, alpha_1, beta_1, current_LL)
      iteration_results[[iteration]] <- result
    }

    return(iteration_results)
  }
  else if (parameterization == "md") {

    # TO DO: code algorithm for mean/dispersion parameterization

  }
  else {

    stop(paste0("Error: valid parameterization not provided for the EM algorithm. Please ensure that ",
                "your parameterization is either 'shape' or 'md' to indicate whether the model fitting ",
                "should be performed using the alpha/beta shape parameterization or mean/dispersion ",
                "parameterization, respectively."))
  }
}

# Beta density/log-likelihood of each mixture model
bm_density_ <- function(Y, alpha_0, beta_0, alpha_1, beta_1, mix_params, log = TRUE) {

  # Log-sum-exp implementation
  log_lower <- log((1 - mix_params)) + dbeta(Y, alpha_0, beta_0, log = TRUE)
  log_upper <- log(mix_params) + dbeta(Y, alpha_1, beta_1, log = TRUE) # matrices of size n x d

  m <- pmax(log_lower, log_upper)
  log_marginal <- m + log(exp(log_lower - m) + exp(log_upper - m))

  if (log) {

    return (log_marginal)
  }
  else {

    return (exp(log_marginal))
  }
}

# E step for shape parametrization
e_step_ <- function(Y, mix_params, alpha_0, beta_0, alpha_1, beta_1) {

  n <- nrow(Y)

  # Log-sum-exp implementation
  log_lower <- log((1 - mix_params)) + dbeta(Y, alpha_0, beta_0, log = TRUE)
  log_upper <- log(mix_params) + dbeta(Y, alpha_1, beta_1, log = TRUE)

  m <- pmax(log_lower, log_upper)
  log_marginal <- m + log(exp(log_lower - m) + exp(log_upper - m))

  upper_responsibility <- exp(log_upper - log_marginal) # Note that lower_responsibility = 1 - upper_responsibility
  mix_params <- 1/n * colSums(upper_responsibility) # Desired vector of length d

  return(mix_params)
}

# M step for shape parameterization
# Reminder that all elements of mix_params are held constant for each run of maxLik (inner loop of EM algorithm)
m_step_ <- function(Y, mix_params, alpha_0, beta_0, alpha_1, beta_1) {

  d <- ncol(Y)

  # d x 4 matrix containing start values for MLE
  start <- cbind(alpha_0, beta_0, alpha_1, beta_1)

  start_list <- lapply(seq_len(d), function(t) {

    v <- as.numeric(start[t, ])
    names(v) <- c("alpha_0", "beta_0", "alpha_1", "beta_1")
    v
  })

  # Likelihood for an individual mixture model
  # Data is called x (not Y) so as to be compatible with maxLik()
  individual_likelihood_ <- function(params, x, mix_param_t) {

    alpha_0 <- params[1] # I personally don't like the numerical indexing, but am too lazy to make params a named numeric right now
    beta_0 <- params[2]
    alpha_1 <- params[3]
    beta_1 <- params[4]

    if(any(c(alpha_0, beta_0, alpha_1, beta_1) <= 0)) {

      return(NA)
    }
    return(bm_density_(x, alpha_0, beta_0, alpha_1, beta_1, mix_param_t, log = TRUE)) # All scalar inputs to bm_density_()
  }

  # Runs maxLik() on an individual mixture model
  individual_fit_ <- function(start_row, Y_col, mix_param_t) {

    alpha_0 <- as.numeric(start_row["alpha_0"])
    beta_0 <- as.numeric(start_row["beta_0"])
    alpha_1 <- as.numeric(start_row["alpha_1"])
    beta_1 <- as.numeric(start_row["beta_1"])

    individual_start <- c(alpha_0, beta_0, alpha_1, beta_1)

    output <- maxLik(logLik = individual_likelihood_,
                     start = individual_start,
                     x = Y_col,
                     mix_param_t = mix_param_t)

    new_params <- output$estimate
    names(new_params) <- c("alpha_0", "beta_0", "alpha_1", "beta_1")
    new_params
  }

  # To apply maxLik() to all models
  updated_shape_params <- mapply(FUN = individual_fit_,
                                 start_row = start_list,
                                 Y_col = as.data.frame(Y),
                                 mix_param_t = mix_params,
                                 SIMPLIFY = FALSE)
  updated_shape_params <- unlist(updated_shape_params)
  updated_shape_params <- t(matrix(updated_shape_params, nrow = 4, ncol = d))

  colnames(updated_shape_params) <- c("alpha_0", "beta_0", "alpha_1", "beta_1")
  rownames(updated_shape_params) <- colnames(Y)

  return(updated_shape_params) # Returns the output formatted as matrix of size d x 4
}
