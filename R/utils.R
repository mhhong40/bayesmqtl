# Beta densities but infinities are replaced with 0's
# dbeta_no_NAs_ <- function(X, alpha, beta) {
#
#   db <- dbeta(X, alpha, beta)
#   db <- ifelse(is.infinite(db), 0, db)
# }
