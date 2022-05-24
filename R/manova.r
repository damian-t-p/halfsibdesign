#' Compute 1- or 2-way MANOVA estimates
#'
#' Given either a tuple of sum-of-squares matrices or a dataframe, computes
#' MANOVA estimates for the covariances
#' @name manova
#' @seealso \code{\link{manova_2way}} for generic version.
NULL
#> NULL

#' @rdname manova
#' @param Mw,Mb Sum-of-squares within and between matrices
#' @param r,n Number of
#' @export
manova_1way <- function(Mw, r, Mb, n) {
  list(
    Sb = (Mb - Mw)/r,
    Sw = Mw
  )
}

#' @rdname manova
#' @param M1,M2,M3 Sum-of-squares matrices within sire, dam and individual groups
#' @param I1,I2,I3 Number of sires, dams per sire and individuals per dam
#' @export
manova_2way_mat <- nest_cov_estimate(manova_1way)

#' @rdname manova
#' @export
#' @param df Data frame for 2-way design
manova_2way_df <- df_cov_estimate(manova_2way_mat)
