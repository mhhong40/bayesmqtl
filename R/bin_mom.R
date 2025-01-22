# This file is part of the `bayesmqtl` package:
#   https://github.com/mhhong40/bayesmqtl
#
#' Generate bins as input to bayesmqtl()
#'
#' @param Y The CpG response data matrix.
#' @param bounds A user-supplied vector of upper bounds for each bin.
#' @return A list containing the indices of each CpG belonging to each bin.
#' @export
generate_bin_ind <- function(Y, bounds) {

  nbins <- length(bounds)
  if(is.unsorted(bounds)) bounds <- order(bounds)

  # Recall that Y is an n x d matrix.
  m_means <- colMeans(Y)
  mean_min <- min(m_means)
  mean_max <- max(m_means)

  if(min(bounds) < mean_min) bounds[1] <- mean_min
  if(max(bounds) > mean_max) bounds[nbins] <- mean_max

  bin_ind <- list()
  bin_ind[[1]] <- which(m_means < bounds[1]) # first bin
  for(k in 2:(nbins-1)) {

    bin_ind[[k]] <- which(m_means >= bounds[k - 1] & m_means < bounds[k])
  }
  bin_ind[[nbins]] <- which(m_means >= bounds[nbins - 1] & m_means <= bounds[nbins]) # last bin
  class(bin_ind) <- "bin_ind"
  bin_ind
}

#' Uses the method of moments to estimate Beta mixture
#' shape parameters given a fixed mixture proportion.
#'
#' @param bin_ind An object of class bin_ind. Used to pool data
#'    from all CpGs belonging in a given bin to estimate their
#'    shared mixture shape parameters.
#' @export
estimate_shape_params <- function(Y, bin_ind, mix_prop) {

  # check length of mix_prop
  nbins <- length(bin_ind)

  alpha_0 <- beta_0 <- alpha_1 <- beta_1 <- c()
  for(k in 1:nbins) {

    pooled <- c(Y[, bin_ind[[k]] ])
    thresh <- quantile(pooled, mix_prop[k])
    lower <- pooled[which(pooled <= thresh)]
    upper <- pooled[which(pooled > thresh)]

    m_l <- mean(lower)
    v_l <- var(lower)

    m_u <- mean(upper)
    v_u <- var(upper)

    # Method of moments estimation
    alpha_0[k] <- m_l * ((m_l * (1 - m_l)) / v_l - 1)
    beta_0[k] <- (1 - m_l) * ((m_l * (1 - m_l)) / v_l - 1)

    alpha_1[k] <- m_u * ((m_u * (1 - m_u)) / v_u - 1)
    beta_1[k] <- (1 - m_u) * ((m_u * (1 - m_u)) / v_u - 1)
  }

  list_shape <- create_named_list_(alpha_0, beta_0, alpha_1, beta_1)
  class(list_shape) <- "list_shape"

  list_shape
}
