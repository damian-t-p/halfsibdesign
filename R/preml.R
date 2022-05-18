#' Compute 1- or 2-way pseudo-REML estimates
#'
#' Given either a tuple of sum-of-squares matrices or a dataframe, computes
#' preudo-REML estimates for the covariances
#' @name preml
NULL
#> NULL

#' @rdname preml
#' @param Mw,Mb Sum-of-squares within and between matrices
#' @param r,n Number of
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
#' @export
preml_2way <- nest_cov_estimate(preml_1way)

#' @rdname preml
#' @export
#' @param df Data frame for 2-way design
preml_2way_df <- df_cov_estimate(preml_2way)
