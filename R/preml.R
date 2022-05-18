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

#' @export
preml_2way <- nest_cov_estimate(preml_1way)

#' @export
preml_2way_df <- df_cov_estimate(preml_2way)
