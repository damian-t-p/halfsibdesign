#' Compute 1- or 2-way pseudo-REML estimates
#'
#' Given either a tuple of sum-of-squares matrices or a dataframe, computes
#' preudo-REML estimates for the covariances.
#'
#'
#' The pseudo-REML estimates are non-negative definite estimates of the relevant
#' covariance matrices as described in \[1\]. These are not the exact REML estimates,
#' but coincide with them in the MANOVA estimates are non-negative definite.
#'
#' @name preml
#'
#' @return A list with entries
#' \itemize{
#' \item in the one-way case: `Sb`, `Sw`, which are the estimates of
#' the between-groups and within-groups covariance matrices.
#' \item in the two-way case: `S1`, `S2` and `S3`, which are the estimates of
#' the sires, dams and individual covariance matrices.
#' }
#'
#' @references \[1\] Y. Amemiyah. "What Should be Done When an Estimated between-Group
#' Covariance Matrix is not Nonnegative Definite?" In *The Americal Statistician*
#' 39.2 (1985).
#'
#' @seealso \code{\link{preml_2way}} for generic version.
NULL
#> NULL

#' @rdname preml
#' @param Mw,Mb Sum-of-squares within and between matrices.
#' @param r,n Number individuals per group and groups.
#'
#' @export
preml_1way <- function(Mw, r, Mb, n) {
  U <- chol(Mw)
  L <- solve(U)

  eigs <- eigen(t(L) %*% Mb %*% L, symmetric = TRUE)
  Lambda_thresh <- diag(pmax(eigs$values - 1, 0), nrow=length(eigs$values)) # (Lambda - I)_+
  Q <- eigs$vectors

  P <- t(U) %*% Q

  Sb <- P %*% Lambda_thresh %*% t(P) / r
  Sw <- ((n-1) * (Mb - r*Sb) + n*(r-1) * Mw) / (n*r - 1)

  list(
    Sb = Sb,
    Sw = Sw
  )
}

#' @rdname preml
#' @param M1,M2,M3 Sum-of-squares matrices within sire, dam and individual groups.
#' @param I1,I2,I3 Number of sires, dams per sire and individuals per dam.
#' @export
preml_2way_mat <- nest_cov_estimate(preml_1way)


#' @rdname preml
#' @export
#' @param df Data frame for 2-way design.
#' @param ... Other arguments.
preml_2way_df <- df_cov_estimate(preml_2way_mat)


