#' Compute the REML criterion
#'
#' Computes the REML criterion for a balanced 2-way design.
#'
#' @param S1,S2,S3 Covariance estimates
#' @param M1,M2,M3 Sum-of-squares matrices
#' @param I1,I2,I3 Numbers of groups
#' @export
reml_crit <- function(M1, S1, I1, M2, S2, I2, M3, S3, I3) {
  p <- nrow(M1)

  n1 <- I3*I2*(I1 - 1)
  n2 <- I3*(I2 - 1)
  n3 <- I3 - 1
  n <- n1 + n2 + n3

  G1 <- M1/n1
  G2 <- M2/n2
  G3 <- M3/n3

  sig1 <- S1
  sig2 <- S1 + I1 * S2
  sig3 <- S1 + I1 * S2 + I1*I2 * S3

  c <- -p*n/2 * log(2*pi) +
    -p/2 * log(I1*I2*I3) +
    #-1/2 * log(det(2*pi*sig3)) +
    -1/2 * (
      n1 * log(det(G1)) +
        n2 * log(det(G2)) +
        n3 * log(det(G3))
    )

  1/2 * (
    n1 * (log(det(solve(sig1) %*% G1)) - sum(diag(solve(sig1) %*% G1))) +
      n2 * (log(det(solve(sig2) %*% G2)) - sum(diag(solve(sig2) %*% G2))) +
      n3 * (log(det(solve(sig3) %*% G3)) - sum(diag(solve(sig3) %*% G3)))
  ) + c

}
