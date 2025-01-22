# This file is part of the `bayesmqtl` package:
#   https://github.com/mhhong40/bayesmqtl
#
# Assorted utility functions for basic preprocessing,
# prevention of over/underflow, and cleanup.
# Adapted from the functions in Hélène's locus>utils.

## Adapted from Hélène's locus>prepare_locus check_annealing_() function.
check_annealing_ <- function(anneal) {

  check_structure_(anneal, "vector", "numeric", 3, null_ok = TRUE)

  if (!is.null(anneal)) {

    check_natural_(anneal[c(1, 3)])
    check_positive_(anneal[2])

    stopifnot(anneal[1] %in% 1:3)

    if (anneal[2] < 1.5)
      stop(paste0("Initial temperature very small. May not be large enough ",
                  "for a successful exploration. Please increase it or select no annealing."))

    if (anneal[3] > 1000)
      stop(paste0("Temperature ladder size very large. This may be unnecessarily ",
                  "computationally demanding. Please decrease it."))

  }

}

check_natural_ <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps | abs(x - round(x)) > eps)) {
    stop(paste0(deparse(substitute(x)),
                " must be natural."))
  }
}

check_positive_ <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps)) {
    err_mess <- paste0(deparse(substitute(x)), " must be positive, greater than ",
                       format(eps, digits = 3), ".")
    if (length(x) > 1) err_mess <- paste0("All entries of ", err_mess)
    stop(err_mess)
  }
}

check_zero_one_ <- function(x){
  if (any(x < 0) | any(x > 1)) {
    err_mess <- paste0(deparse(substitute(x)), " must lie between 0 and 1.")
    if (length(x) > 1) err_mess <- paste0("All entries of ", err_mess)
    stop(err_mess)
  }
}

check_structure_ <- function(x, struct, type, size = NULL,
                             null_ok = FALSE,  inf_ok = FALSE, na_ok = FALSE) {
  if (type == "double") {
    bool_type <-  is.double(x)
    type_mess <- "double-precision "
  } else if (type == "integer") {
    bool_type <- is.integer(x)
    type_mess <- "integer "
  } else if (type == "numeric") {
    bool_type <- is.numeric(x)
    type_mess <- "numeric "
  } else if (type == "logical") {
    bool_type <- is.logical(x)
    type_mess <- "boolean "
  } else if (type == "string") {
    bool_type <- is.character(x)
    type_mess <- "string "
  }

  bool_size <- TRUE # for case size = NULL (no assertion on the size/dimension)
  size_mess <- ""
  if (struct == "vector") {
    bool_struct <- is.vector(x) & (length(x) > 0) # not an empty vector
    if (!is.null(size)) {
      bool_size <- length(x) %in% size
      size_mess <- paste0(" of length ", paste0(size, collapse=" or "))
    }
  } else if (struct == "matrix") {
    bool_struct <- is.matrix(x) & (length(x) > 0) # not an empty matrix
    if (!is.null(size)) {
      bool_size <- all(dim(x) == size)
      size_mess <- paste0(" of dimension ", size[1], " x ", size[2])
    }
  }

  correct_obj <- bool_struct & bool_type & bool_size

  bool_null <- is.null(x)

  if (!is.list(x) & type != "string") {
    na_mess <- ""
    if (!na_ok) {
      if (!bool_null) correct_obj <- correct_obj & !any(is.na(x))
      na_mess <- " without missing value"
    }

    inf_mess <- ""
    if (!inf_ok) {
      if (!bool_null) correct_obj <- correct_obj & all(is.finite(x[!is.na(x)]))
      inf_mess <- ", finite"
    }
  } else {
    na_mess <- ""
    inf_mess <- ""
  }

  null_mess <- ""
  if (null_ok) {
    correct_obj <- correct_obj | bool_null
    null_mess <- " or must be NULL"
  }

  if(!(correct_obj)) {
    stop(paste0(deparse(substitute(x)), " must be a non-empty ", type_mess, struct,
                size_mess, inf_mess, na_mess, null_mess, "."))
  }
}

create_named_list_ <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}

get_annealing_ladder_ <- function(anneal, verbose) {

  # ladder set following:
  # Importance Tempering, Robert B. Gramacy & Richard J. Samworth, pp.9-10, arxiv v4

  k_m <- 1 / anneal[2]
  m <- anneal[3]

  if(anneal[1] == 1) {

    type <- "geometric"

    delta_k <- k_m^(1 / (1 - m)) - 1

    ladder <- (1 + delta_k)^(1 - m:1)

  } else if (anneal[1] == 2) { # harmonic spacing

    type <- "harmonic"

    delta_k <- ( 1 / k_m - 1) / (m - 1)

    ladder <- 1 / (1 + delta_k * (m:1 - 1))

  } else { # linear spacing

    type <- "linear"

    delta_k <- (1 - k_m) / (m - 1)

    ladder <- k_m + delta_k * (1:m - 1)
  }

  if (verbose)
    cat(paste0("** Annealing with ", type," spacing ** \n\n"))

  ladder

}

# Directly adapted from Hélène's code for locus>utils>inv_mills_ratio_matrix_
inv_mills_ratio_matrix_ <- function(U, Y) {

  if (is.matrix(U)) m <- matrix(NA, nrow = nrow(U), ncol = ncol(U))
  else m <- rep(NA, length(U))

  U_1 <- U[Y==1]
  m_1 <- exp(dnorm(U_1, log = TRUE) - pnorm(U_1, log.p = TRUE))
  m_1[m_1 < -U_1] <- -U_1

  m[Y==1] <- m_1


  U_0 <- U[Y==0]
  m_0 <- - exp(dnorm(U_0, log = TRUE) - pnorm(U_0, lower.tail = FALSE, log.p = TRUE))
  m_0[m_0 > -U_0] <- -U_0

  m[Y==0] <- m_0

  m

}

#' @export
estimate_mem_probs_ <- function(Y, z, mix_prop, alpha_0, beta_0, alpha_1, beta_1) {

  d <- ncol(Y)
  eps <- .Machine$double.eps^0.5 # to prevent under/overflow when the membership probability equals 0 or 1
                                 # (though, realistically, this should never happen)
  for (t in 1:d) {

    z[, t] <- mix_prop[t]*dbeta(Y[, t], shape1 = alpha_1[t], shape2 = beta_1[t])/((1 - mix_prop[t])*dbeta(Y[, t], shape1 = alpha_0[t], shape2 = beta_0[t]) +
                                                           mix_prop[t]*dbeta(Y[, t], shape1 = alpha_1[t], shape2 = beta_1[t]))
  }

  z[z==0] <- eps
  z[z==1] <- 1 - eps

  z
}

log_one_plus_exp_ <- function(x) { # computes log(1 + exp(x)) avoiding
                                   # numerical overflow
  m <- x
  m[x < 0] <- 0

  log(exp(x - m) + exp(- m)) + m
}

checkpoint_ <- function(it, checkpoint_path, theta_0_vb, theta_1_vb,
                        converged, lb_new, lb_old, rate = 100) {

  if (!is.null(checkpoint_path) && it %% rate == 0) {

    diff_lb <- abs(lb_new - lb_old)

    tmp_vb <- create_named_list_(theta_0_vb, theta_1_vb, converged, it, lb_new, diff_lb)

    file_save <- paste0(checkpoint_path, "tmp_output_it_", it, ".RData")

    save(tmp_vb, file = file_save)

    old_file_clean_up <- paste0(checkpoint_path, "tmp_output_it_", it - 2 * rate, ".RData") # keep only the last two for comparison

    if (file.exists(old_file_clean_up))
      file.remove(old_file_clean_up)

  }

}

checkpoint_clean_up_ <- function(checkpoint_path) {

  if (!is.null(checkpoint_path)) {

    old_files_clean_up <- list.files(path = checkpoint_path, pattern = "tmp_output_it_")

    sapply(old_files_clean_up, function(ff) {
      if (file.exists(file.path(checkpoint_path, ff)))
        file.remove(file.path(checkpoint_path, ff))
    })

  }

}
