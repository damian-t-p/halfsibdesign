#' Simulate a dataset from a full-sib design
#'
#' @param mu Mean vector.
#' @param A,E Covariance matrices of the sires dams and individuals.
#' @param I,J Number of sires, dams per sire and individuals per dam.
#'
#' @export
rfullsib <- function(mu, A, I, E, J) {
  q <- length(mu)

  traits <- rep(1:q, times = I*J)
  mu_full = rep(mu, times = I*J)

  zero_mean <- rep(0, q)

  sires <- rep(1:I, each = J*q)
  alphas <- MASS::mvrnorm(I, zero_mean, A) %>%
    hrepmat(J) %>%
    {c(t(.))}

  individuals <- rep(1:(I*J), each = q)
  epsilons <- MASS::mvrnorm(I*J, zero_mean, E) %>%
    {c(t(.))}

  tibble::tibble(
    trait = factor(traits),
    value = mu_full + alphas + epsilons,
    sire = factor(sires),
    individual = factor(individuals)
  )
}
