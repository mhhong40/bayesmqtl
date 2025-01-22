# This file is part of the `bayesmqtl` package:
#   https://github.com/mhhong40/bayesmqtl
#
# Contains code for basic data preparation functions and sanity checks.
# Much of it is adapted from Hélène's locus>prepare_locus functions.

## Adapted from Hélène's locus>prepare_locus prepare_blocks() and
#  internal set_blocks() functions. Used for splitting the data into
#  predictor blocks for effective parallelization.
prepare_blocks <- function(list_blocks) {

  if (!inherits(list_blocks, "blocks"))
    stop(paste0("The provided list_blocks must be an object of class ``blocks''. \n",
                "*** you must either use the function set_blocks to give the settings ",
                "for parallel applications of bayesmqtl on blocks of candidate ",
                "effector SNPs or set list_blocks to NULL to apply locus jointly on ",
                "all the candidate effectors (sufficient RAM required). ***"))

  vec_fac_bl <- list_blocks$vec_fac_bl

  tab_bl <- table(vec_fac_bl)
  pres_bl <- tab_bl > 0

  n_bl  <- sum(pres_bl)
  if(list_blocks$n_cpus > n_bl) n_cpus <- n_bl
  else n_cpus <- list_blocks$n_cpus

  create_named_list_(n_bl, n_cpus, vec_fac_bl)
}

#' @export
set_blocks <- function(p, pos_bl, n_cpus, verbose = TRUE) {

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p) # from utils

  check_structure_(verbose, "vector", "logical", 1)

  check_structure_(pos_bl, "vector", "numeric")
  check_natural_(pos_bl)

  if (length(pos_bl) > 25)
    warning(paste0("The provided number of blocks may be too large for accurate ",
                   "inference. If possible, use less blocks."))

  if (any(pos_bl < 1) | any(pos_bl > p))
    stop("The positions provided in pos_bl must range between 1 and total number of variables in X, p.")

  if (any(duplicated(pos_bl)))
    stop("The positions provided in pos_bl must be unique.")

  if (any(pos_bl != cummax(pos_bl)))
    stop("The positions provided in pos_bl must be monotonically increasing.")

  vec_fac_bl <- as.factor(cumsum(seq_along(1:p) %in% pos_bl))

  n_bl <- length(unique(vec_fac_bl))

  check_structure_(n_cpus, "vector", "numeric", 1)
  check_natural_(n_cpus)

  if (n_cpus > 1) {

    n_cpus_avail <- parallel::detectCores()
    if (n_cpus > n_cpus_avail) {
      n_cpus <- n_cpus_avail
      warning(paste0("The number of CPUs specified exceeds the number of CPUs ",
                     "available on the machine. The latter has been used instead."))
    }

    if (n_cpus > n_bl){
      message <- paste0("The number of cpus in use is at most equal to the number of blocks.",
                        "n_cpus is therefore set to ", n_bl, ". \n")
      if(verbose) cat(message)
      else warning(message)
      n_cpus <- n_bl
    }

    if (verbose) cat(paste0("bayesmqtl applied in parallel on ", n_bl,
                            " blocks of candidate effectors, using ", n_cpus, " CPUs.\n",
                            "Please make sure that enough RAM is available. \n"))
  }

  p_blocks <- p

  list_blocks <- create_named_list_(p_blocks, n_bl, n_cpus, vec_fac_bl)

  class(list_blocks) <- "blocks"

  list_blocks
}

prepare_list_hyper_ <- function(list_hyper, Y, p, names_x, names_y, verbose) {

  d <- ncol(Y)
  p_hyper_match <- p

  if (is.null(list_hyper)) {

    if (verbose) stop(paste0("There is no default hyperparameter setting function at the moment.",
                     " Please specify a hyperparameter list."))

  }
  else {

    if (!inherits(list_hyper, c("hyper", "out_hyper")))
      stop(paste0("The provided list_hyper must be an object of class ``hyper'' ",
                  "or ``out_hyper''. \n",
                  "*** you must either use the function set_hyper to ",
                  "set your own hyperparameters or use list_hyper from a ``vb`` ",
                  "object. ***"))

    if (list_hyper$d_hyper != d)
      stop(paste0("The dimensions (d) of the provided hyperparameters ",
                  "(list_hyper) are not consistent with that of Y.\n"))

    if (list_hyper$p_hyper != p_hyper_match)
      stop(paste0("The dimensions (p) of the provided hyperparameters ",
                  "(list_hyper) are not consistent with that of X.\n"))

  }

  class(list_hyper) <- "out_hyper"

  list_hyper
}

prepare_list_init_ <- function(list_init, Y, p, user_seed, verbose) {

  d <- ncol(Y)
  p_init_match <- p

  if (is.null(list_init)) {

    if (verbose) stop(paste0("There is no default parameter setting function at the moment.",
                             " Please specify a parameter list."))

  }
  else {

    if (!is.null(user_seed))
      warning("user_seed not used since a non-NULL list_init was provided. \n")

    if (!inherits(list_init, c("init", "out_init")))
      stop(paste0("The provided list_init must be an object of class ``init`` or ",
                  " `` out_init``. \n",
                  "*** you must either use the function set_init to ",
                  "set your own initialization or use list_init from a ``vb`` ",
                  "object. ***"))

    if (list_init$d_init != d)
      stop(paste0("The dimensions (d) of the provided initial parameters ",
                  "(list_init) are not consistent with that of Y.\n"))

    if (list_init$p_init != p_init_match)
      stop(paste0("The dimensions (p) of the provided initial parameters ",
                  "(list_init) are not consistent with that of X.\n"))

  }
  class(list_init) <- "out_init"

  list_init
}

prepare_bin_ind_ <- function(bin_ind, Y, verbose) {

  d <- ncol(Y)
  k <- length(bin_ind)

  if (is.null(bin_ind)) {

    if (verbose) stop(paste0("A pre-generated bin index list is required for the ",
                             "estimation of the remaining model parameters. \n",
                             "Please specify a bin index list."))

  }
  else {
    if (!inherits(bin_ind, c("bin_ind")))
      stop(paste0("The provided bin_ind must be an object of class ``bin_ind``. \n",
                  "*** You must use the function generate_bin_ind to obtain ",
                  "bin indexes for each response variable to use as \n",
                  "input to properly fetch their mixture components' shape parameters. ***"))

    total <- sum(unlist(lapply(bin_ind, sum)))
    expected <- d * (d + 1) / 2
    if(total > expected)
      stop(paste0("Each CpG must be assigned to exactly one bin. \n.",
                  "Please check for duplicate indices in the bin_ind object."))
    else if(total < expected)
      stop(paste0("When combined, the bins must cover every CpG. \n",
                  "Please check for missing indices in the bin_ind object using ",
                  "which(1:d %in% unlist(bin_ind) == FALSE)."))
  }
}

prepare_list_shape_ <- function(list_shape, Y, k, verbose) {

  d <- ncol(Y)

  if (is.null(list_shape)) {

    if (verbose) stop(paste0("A pre-generated shape parameter list is required for the ",
                             "estimation of the remaining model parameters. \n",
                             "Please specify a parameter list."))

  }
  else {
    if (!inherits(list_shape, c("list_shape")))
      stop(paste0("The provided list_shape must be an object of class ``list_shape``. \n",
                  "*** You must use the function estimate_shape_params to obtain ",
                  "bin-wide shape parameter estimates that can be used as \n",
                  "input to the VI algorithm for the remaining parameters. ***"))
    if (!all(lengths(list_shape) == k)) {

      stop(paste0("Each shape parameter vector must be of length k = ", k,
                  " for proper assignment. Please ensure that this condition is met."))
    }
  }
}

prepare_data_ <- function(Y, X, user_seed, tol, maxit, verbose, checkpoint_path) {

  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)

  check_structure_(tol, "vector", "numeric", 1)
  check_positive_(tol, eps=.Machine$double.eps)

  check_structure_(maxit, "vector", "numeric", 1)
  check_natural_(maxit)

  check_structure_(X, "matrix", "numeric")

  if (!is.null(checkpoint_path)) {

    if (!dir.exists(checkpoint_path)) {
      stop(paste0("The directory specified in checkpoint_path doesn't exist. ",
                  "Please make sure to provide a valid path."))
    }
  }

  n <- nrow(X)
  p <- ncol(X)

  check_structure_(Y, "matrix", "numeric")
  d <- ncol(Y)

  if (nrow(Y) != n) stop("X and Y must have the same number of samples.")

  if (is.null(rownames(X)) & is.null(rownames(Y)))
    rownames(X) <- rownames(Y) <- paste0("Ind_", 1:n)
  else if (is.null(rownames(X))) rownames(X) <- rownames(Y)
  else if (is.null(rownames(Y))) rownames(Y) <- rownames(X)
  else if (any(rownames(X) != rownames(Y)))
    stop("The provided rownames of X and Y must be the same.")

  if (is.null(colnames(X))) colnames(X) <- paste0("Eff_x_", 1:p)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("Resp_", 1:d)

  p <- ncol(X)

  create_named_list_(Y, X)

}
