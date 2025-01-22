# function tests...

# updateShapeA in R
update_a_ <- function(a_init, a, b, c, d, c_temp = 1) {

  alpha_bar <- a/b
  beta_bar <- c/d
  a <- c_temp * (a_init + alpha_bar) * (digamma(alpha_bar + beta_bar) - digamma(alpha_bar)
                                        + beta_bar * trigamma(alpha_bar + beta_bar) *
                                          ( (digamma(c) - log(d)) - log(beta_bar)) )
}

update_c_ <- function(c_init, a, b, c, d, c_temp = 1) {

  alpha_bar <- a/b
  beta_bar <- c/d
  c = c_temp * (c_init + beta_bar) * (digamma(alpha_bar + beta_bar) - digamma(beta_bar)
                                               + alpha_bar * trigamma(alpha_bar + beta_bar) *
                                                 ( (digamma(a) - log(b)) - log(alpha_bar)) )
}

computeElogDiffSq <- function(x) {

  (digamma(x) - log(x))^2 + trigamma(x)
}

compute_p_ <- function(a, b, c, d, alpha_k, beta_k) {

  alpha_bar = a/b
  beta_bar = c/d

  Elogalpha = digamma(a) - log(b)
  Elogbeta = digamma(c) - log(d)

  P = log((1.0 / beta(alpha_k, beta_k))) +
    alpha_bar * (digamma(alpha_k + beta_k) - digamma(alpha_bar)) * (Elogalpha - log(alpha_bar)) +
    beta_bar * (digamma(alpha_k + beta_k) - digamma(beta_bar)) * (Elogbeta - log(beta_bar)) +
    0.5 * alpha_bar^2 * (trigamma(alpha_bar + beta_bar) - trigamma(alpha_bar)) * computeElogDiffSq(a) +
    0.5 * beta_bar^2 * (trigamma(alpha_bar + beta_bar) - trigamma(beta_bar)) * computeElogDiffSq(c) +
    alpha_bar * beta_bar * trigamma(alpha_bar + beta_bar) * (Elogalpha - log(alpha_bar)) * (Elogbeta - log(beta_bar));
}
