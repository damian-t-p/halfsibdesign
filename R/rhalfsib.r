#' Simulate a dataset from a half-sib design
#'
#' Sample responses from the following model:
#' \deqn{y_{ijk} = \mu + \alpha_i + \beta_{ij} + \epsilon_{ijk},}
#' where \eqn{\mu} is a deterministic mean vector and each \eqn{\alpha_i,
#' \beta_{ij}, \epsilon_{ijk}} are independent mean-zero multivariate
#' Gaussians with covariance matrices \eqn{A, B, E} respectively.
#'
#' @param mu Mean vector.
#' @param A,B,E Covariance matrices of the sires dams and individuals.
#' @param I,J,K Number of sires, dams per sire and individuals per dam.
#'
#' @return A `tibble` with column names `trait`, `value`, `sire`,
#' `dam`, `individual`. The `value` column contains the response values whose
#' dimension is indexed by `trait`.
#'
#' @export
rhalfsib <- function(mu, A, I, B, J, E, K) {
  q <- length(mu)

  traits <- rep(1:q, times = I*J*K)
  mu_full = rep(mu, times = I*J*K)

  zero_mean <- rep(0, q)

  sires <- rep(1:I, each = J*K*q)
  alphas <- MASS::mvrnorm(I, zero_mean, A) %>%
    hrepmat(J*K) %>%
    {c(t(.))}

  dams <- rep(1:(I*J), each = K*q)
  betas <- MASS::mvrnorm(I*J, zero_mean, B) %>%
    hrepmat(K) %>%
    {c(t(.))}

  individuals <- rep(1:(I*J*K), each = q)
  epsilons <- MASS::mvrnorm(I*J*K, zero_mean, E) %>%
    {c(t(.))}

  tibble::tibble(
    trait = factor(traits),
    value = mu_full + alphas + betas + epsilons,
    sire = factor(sires),
    dam = factor(dams),
    individual = factor(individuals)
  )
}

hrepmat <- function(X, n) {
  matrix(rep(X, n), nrow = nrow(X))
}
