#' Contains the various functions needed to fit two-component Beta mixture models
#' to DNA methylation (microarray) data.
#'
#' Fits the mixture model via an expectation-maximization (EM) algorithm.
#' E and M-steps are both coded from scratch for numerical stability (using log-sum-exp) and to facilitate parallel execution.
#' @param Y The methylation data matrix of size n x d, where n = the number of samples and d = the number of methylation sites.
#' @param threshold A value ranging from 0 to 1 (exclusive) indicating how the data distribution should first be split for parameter initialization. Default is 0.5 (data is split into upper and lower halves).
#' @param parametrization A string that is either "shape" or "md." Specifying "shape" performs the model fitting based on the alpha/beta shape parametrization of the Beta distribution, while md performs that based on the mean/dispersion parametrization.
#' @param tol Tolerance value for the difference in the previous and current iterations' log-likelihood values that determines when the EM algorithm should stop. Default is 0.1.
#' @param maxit Maximum allowed number of iterations for the EM algorithm before it force quits.
#' @param obj_param Indicates which of the five mixture model parameters is sought to be estimated. If supplied with "all," then the true_params argument should be supplied with NULL (the EM algorithm estimates all parameters entirely from scratch).
#' @param true_params A list of values supplying the true values of all other model parameters aside from that supplied to obj_param. List elements must have the model parameter names; for instance, if obj_param == "pi", then true_params contains alpha_0, beta_0, etc.
#' @param sens_param (Assumes that obj_param == "all.") If not "none," specifies the shape parameter for which various initializations are being supplied for sensitivity analysis. If "none," the default method-of-moments initialization is used for all shape parameters.
#' @param sens_param_val A vector containing initial values for the shape parameter specified by the sens_param argument.
fit_mixture_model_ <- function(Y, threshold = 0.5, parametrization = c("shape", "md"), tol = 0.1, maxit = 1000,
                               obj_param = c("all", "pi", "alpha_0", "beta_0", "alpha_1", "beta_1"),
                               sens_param = c("none", "alpha_0", "beta_0", "alpha_1", "beta_1"),
                               true_params, sens_param_val) {

  n <- nrow(Y)
  d <- ncol(Y)

  # Control over numerical under/overflow
  eps <- 1e-6
  Y <- pmin(pmax(Y, eps), 1 - eps)

  cutoffs <- matrixStats::colQuantiles(Y, probs = threshold)

  lower_mask <- sweep(Y, 2, cutoffs, '>=')
  lower <- Y * (1 - lower_mask)
  upper <- Y * lower_mask

  # Initialize model parameters (method of moments initialization)
  m_l <- matrixStats::colMeans2(lower)
  v_l <- matrixStats::colVars(lower)

  m_u <- matrixStats::colMeans2(upper)
  v_u <- matrixStats::colVars(upper)

  if (parametrization == "shape") {

    # Initialize model parameters (method of moments initialization for shape parameters)
    if(obj_param == "all") {

      pi <- rep(threshold, d)

      # Check if any shape parameters have initializations directly supplied to them for sensitivity analysis
      if(sens_param == "none") {

        alpha_0 <- m_l * (m_l * (1 - m_l) / v_l - 1)
        beta_0 <- (1 - m_l) * (m_l * (1 - m_l) / v_l - 1)
        alpha_1 <- m_u * (m_u * (1 - m_u) / v_u - 1)
        beta_1 <- (1 - m_u) * (m_u * (1 - m_u) / v_u - 1)
      }
      else if(sens_param == "alpha_0") {

        alpha_0 <- sens_param_val
        beta_0 <- (1 - m_l) * (m_l * (1 - m_l) / v_l - 1)
        alpha_1 <- m_u * (m_u * (1 - m_u) / v_u - 1)
        beta_1 <- (1 - m_u) * (m_u * (1 - m_u) / v_u - 1)
      }
      else if(sens_param == "beta_0") {

        alpha_0 <- m_l * (m_l * (1 - m_l) / v_l - 1)
        beta_0 <- sens_param_val
        alpha_1 <- m_u * (m_u * (1 - m_u) / v_u - 1)
        beta_1 <- (1 - m_u) * (m_u * (1 - m_u) / v_u - 1)
      }
      else if(sens_param == "alpha_1") {

        alpha_0 <- m_l * (m_l * (1 - m_l) / v_l - 1)
        beta_0 <- (1 - m_l) * (m_l * (1 - m_l) / v_l - 1)
        alpha_1 <- sens_param_val
        beta_1 <- (1 - m_u) * (m_u * (1 - m_u) / v_u - 1)
      }
      else {

        alpha_0 <- m_l * (m_l * (1 - m_l) / v_l - 1)
        beta_0 <- (1 - m_l) * (m_l * (1 - m_l) / v_l - 1)
        alpha_1 <- m_u * (m_u * (1 - m_u) / v_u - 1)
        beta_1 <- sens_param_val
      }
    }
    else if(obj_param == "pi") {

      pi <- rep(threshold, d)
      alpha_0 <- true_params$alpha_0
      beta_0 <- true_params$beta_0
      alpha_1 <- true_params$alpha_1
      beta_1 <- true_params$beta_1
    }
    else if(obj_param == "alpha_0") {

      pi <- true_params$pi
      alpha_0 <- m_l * (m_l * (1 - m_l) / v_l - 1)
      beta_0 <- true_params$beta_0
      alpha_1 <- true_params$alpha_1
      beta_1 <- true_params$beta_1
    }
    else if(obj_param == "beta_0") {

      pi <- true_params$pi
      alpha_0 <- true_params$alpha_0
      beta_0 <- (1 - m_l) * (m_l * (1 - m_l) / v_l - 1)
      alpha_1 <- true_params$alpha_1
      beta_1 <- true_params$beta_1
    }
    else if(obj_param == "alpha_1") {

      pi <- true_params$pi
      alpha_0 <- true_params$alpha_0
      beta_0 <- true_params$beta_0
      alpha_1 <- m_u * (m_u * (1 - m_u) / v_u - 1)
      beta_1 <- true_params$beta_1
    }
    else {

      pi <- true_params$pi
      alpha_0 <- true_params$alpha_0
      beta_0 <- true_params$beta_0
      alpha_1 <- true_params$alpha_1
      beta_1 <- (1 - m_u) * (m_u * (1 - m_u) / v_u - 1)
    }

    current_LL <- colSums(bm_density_(Y, alpha_0, beta_0, alpha_1, beta_1, pi, log = TRUE))
    prev_LL <- current_LL - 100 # arbitrary starting point

    result <- cbind(pi, alpha_0, beta_0, alpha_1, beta_1, current_LL)

    iteration <- 1
    iteration_results <- list()
    iteration_results[[iteration]] <- result

    # Only run EM algorithm on CpGs that have not yet converged
    not_converged <- rep(TRUE, d)

    while(length(not_converged) > 0) { # Runs the algorithm until all models have reached convergence

      if(any(current_LL < prev_LL)) { # For debugging

        stop(paste0("Warning: at least one of the model log-likelihoods is not monotonically increasing with each iteration. Check your implementation."))
      }

      if(iteration == maxit) {

        paste0("Warning: the maximum allowed number of iterations (", maxit, ") was reached before the EM algorithm reached convergence on all CpGs. ",
               "Consider first filtering out particularly low-quality datasets and/or choosing a greater tolerance.")

        return(iteration_results)
      }
      iteration <- iteration + 1
      start_time <- proc.time()[[3]]

      # E-step is always the same
      resp <- e_step_(Y = Y[, not_converged, drop = FALSE],
                      pi = pi[not_converged],
                      alpha_0 = alpha_0[not_converged], beta_0 = beta_0[not_converged],
                      alpha_1 = alpha_1[not_converged], beta_1 = beta_1[not_converged])

      if(obj_param == "pi") {

        # Only pi is optimized via the M-step
        updated_params <- m_step_(Y = Y[, not_converged, drop = FALSE], resp = resp[, not_converged, drop = FALSE],
                                  alpha_0 = alpha_0[not_converged], beta_0 = beta_0[not_converged],
                                  alpha_1 = alpha_1[not_converged], beta_1 = beta_1[not_converged],
                                  obj_param = obj_param)
      }
      else if(obj_param != "all") {

        # The resp argument isn't used since we know the true values of pi
        updated_params <- m_step_(Y[, not_converged, drop = FALSE], pi = pi[not_converged],
                                  alpha_0 = alpha_0[not_converged], beta_0 = beta_0[not_converged],
                                  alpha_1 = alpha_1[not_converged], beta_1 = beta_1[not_converged],
                                  obj_param = obj_param)
      }
      else {

       updated_params <- m_step_(Y[, not_converged, drop = FALSE], resp = resp[, not_converged, drop = FALSE],
                                 alpha_0 = alpha_0[not_converged], beta_0 = beta_0[not_converged],
                                 alpha_1 = alpha_1[not_converged], beta_1 = beta_1[not_converged],
                                 obj_param = obj_param)
      }

      pi[not_converged] <- updated_params[, 1]
      alpha_0[not_converged] <- updated_params[, 2]
      beta_0[not_converged] <- updated_params[, 3]
      alpha_1[not_converged] <- updated_params[, 4]
      beta_1[not_converged] <- updated_params[, 5]

      prev_LL <- current_LL
      current_LL <- colSums(bm_density_(Y, alpha_0, beta_0, alpha_1, beta_1, pi, log = TRUE)) # With updated parameters

      end_time <- proc.time()[[3]]
      time_elapsed <- end_time - start_time

      # Update the indices of CpGs that still have not converged
      not_converged <- abs(current_LL - prev_LL) > tol
      not_converged <- which(not_converged)

      # Storing step results and runtime
      updated_params <- cbind(updated_params, current_LL)
      iteration_results[[iteration]] <- list(updated_params, time_elapsed)
      names(iteration_results[[iteration]]) <- c("param_estimates", "time_elapsed") # To access runtimes only, use unlist(sapply(fit_results, function(x) x[2]))
    }
    return(iteration_results)
  }
  else if (parametrization == "md") {

    # TO DO: code algorithm for mean/dispersion parametrization

  }
  else {

     stop(paste0("Error: valid parametrization not provided for the EM algorithm. Please ensure that ",
                "your parametrization is either 'shape' or 'md' to indicate whether the model fitting ",
                "should be performed using the alpha/beta shape parametrization or mean/dispersion ",
                "parametrization, respectively."))
  }
}

# E step for shape parametrization
e_step_ <- function(Y, pi, alpha_0, beta_0, alpha_1, beta_1) {

  n <- nrow(Y)

  # Log-sum-exp implementation
  log_lower <- sweep(dbeta(Y, alpha_0, beta_0, log = TRUE), 2, log1p(-pi), "+") # log1p(-pi) for added numerical stability
  log_upper <- sweep(dbeta(Y, alpha_1, beta_1, log = TRUE), 2, log(pi), "+")

  m <- pmax(log_lower, log_upper)
  log_marginal <- m + log(exp(log_lower - m) + exp(log_upper - m))

  upper_responsibility <- exp(log_upper - log_marginal) # Note that lower_responsibility = 1 - upper_responsibility

  return(upper_responsibility) # This is a matrix of size n x d
}

# Beta density/log-likelihood of each mixture model
bm_density_ <- function(Y, alpha_0, beta_0, alpha_1, beta_1, pi, log = TRUE) {

  # Log-sum-exp implementation
  log_lower <- sweep(dbeta(Y, alpha_0, beta_0, log = TRUE), 2, log1p(-pi), "+") # log1p(-pi) for added numerical stability
  log_upper <- sweep(dbeta(Y, alpha_1, beta_1, log = TRUE), 2, log(pi), "+")

  m <- pmax(log_lower, log_upper)
  log_marginal <- m + log(exp(log_lower - m) + exp(log_upper - m))

  if (log) {

    return(log_marginal)
  }
  else {

    return (exp(log_marginal))
  }
}

# Operates on one CpG and is used as input to individual_likelihood_(), hence the data being passed as x
bm_LL_ <- function(x, alpha_0, beta_0, alpha_1, beta_1, pi_t) {

  # Log-sum-exp implementation
  log_lower <- log1p(-pi_t) + dbeta(x, alpha_0, beta_0, log = TRUE)
  log_upper <- log(pi_t) + dbeta(x, alpha_1, beta_1, log = TRUE)

  m <- pmax(log_lower, log_upper)
  log_marginal <- m + log(exp(log_lower - m) + exp(log_upper - m))

  return(sum(log_marginal))
}

# M step for shape parametrization
# When obj_param is specified to be one of the shape parameters, only obj_param is overwritten by MLE; all others are fixed
# The resp argument is only used when either pi or all parameters are the objective
# The pi argument is only used when one of the shape parameters is the objective
m_step_ <- function(Y, resp, pi, alpha_0, beta_0, alpha_1, beta_1, obj_param) {

  eps <- 1e-6 # Can be tweaked, but a Beta mixture is likely extremely misspecified if any parameter reaches such a low value
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
  individual_likelihood_ <- function(params, x, pi_t) {

    alpha_0 <- params[1] # I personally don't like the numerical indexing, but am too lazy to make params a named numeric right now
    beta_0 <- params[2]
    alpha_1 <- params[3]
    beta_1 <- params[4]

    if(any(c(alpha_0, beta_0, alpha_1, beta_1) <= eps)) {

      return(NA)
    }
    return(bm_LL_(x, alpha_0, beta_0, alpha_1, beta_1, pi_t))
    # return(sum(log((1 - pi_t)*dbeta(x, alpha_0, beta_0) + pi_t*dbeta(x, alpha_1, beta_1)))) <- potentially unstable
  }

  if(obj_param == "pi" | obj_param == "all") {

    pi <- matrixStats::colMeans2(resp) # Desired vector of length d; using the matrixStats implementation of this because it is fully optimized

    if(obj_param == "pi") {

      new_params <- cbind(pi, start)
      return(new_params) # Exit the M-step immediately since the shape parameters are fixed; when obj_param is anything else, we proceed to the remainder of the M-step code
    }
  }

  # Runs maxLik() on an individual mixture model when the objective parameter is not pi
  individual_fit_ <- function(start_row, Y_col, pi_t, obj_param = obj_param) {

    alpha_0 <- as.numeric(start_row["alpha_0"])
    beta_0 <- as.numeric(start_row["beta_0"])
    alpha_1 <- as.numeric(start_row["alpha_1"])
    beta_1 <- as.numeric(start_row["beta_1"])

    individual_start <- c(alpha_0, beta_0, alpha_1, beta_1)
    names(individual_start) <- c("alpha_0", "beta_0", "alpha_1", "beta_1")

    if(obj_param == "all") {

      output <- maxLik(logLik = individual_likelihood_,
                       start = individual_start,
                       x = Y_col,
                       pi_t = pi_t) # None of the shape parameters are fixed
    }
    else if(obj_param == "alpha_0") {

      fixed <- c("beta_0", "alpha_1", "beta_1")
    }
    else if(obj_param == "beta_0") {

      fixed <- c("alpha_0", "alpha_1", "beta_1")
    }
    else if(obj_param == "alpha_1") {

      fixed <- c("alpha_0", "beta_0", "beta_1")
    }
    else if(obj_param == "beta_1") {

      fixed <- c("alpha_0", "beta_0", "alpha_1")
    }

    if(exists("fixed")) {

      output <- maxLik(logLik = individual_likelihood_,
                       start = individual_start,
                       x = Y_col,
                       pi_t = pi_t, fixed = fixed)
    }

    new_params <- output$estimate
    names(new_params) <- c("alpha_0", "beta_0", "alpha_1", "beta_1")
    new_params
  }

  # To apply maxLik() to all models
  updated_shape_params <- mapply(FUN = individual_fit_,
                                 start_row = start_list,
                                 Y_col = as.data.frame(Y),
                                 pi_t = pi,
                                 obj_param = rep(obj_param, d),
                                 SIMPLIFY = FALSE)
  updated_shape_params <- unlist(updated_shape_params)
  updated_shape_params <- t(matrix(updated_shape_params, nrow = 4, ncol = d))

  # TO DO: Enforce labels
  updated_shape_params <- enforce_labels_(updated_shape_params) # function is in utils.R
  updated_params <- cbind(pi, updated_shape_params)

  colnames(updated_params) <- c("pi", "alpha_0", "beta_0", "alpha_1", "beta_1")
  rownames(updated_params) <- colnames(Y)

  return(updated_params) # Returns the output formatted as matrix of size d x 4
}
