#' Compute 1- or 2-way MANOVA estimates
#'
#' Given either a tuple of sum-of-squares matrices or a dataframe, computes
#' MANOVA estimates for the covariances
#'
#' @name manova
#'
#' @return A list with entries
#' \itemize{
#' \item in the one-way case: `Sb`, `Sw`, which are the estimates of
#' the between-groups and within-groups covariance matrices.
#' \item in the two-way case: `S1`, `S2` and `S3`, which are the estimates of
#' the sires, dams and individual covariance matrices.
#' }
#'
#' @seealso \code{\link{manova_2way}} for generic version.
NULL
#> NULL

#' @rdname manova
#' @param Mw,Mb Sum-of-squares within and between matrices.
#' @param r,n Number individuals per group and groups.
#' @export
manova_1way <- function(Mw, r, Mb, n) {
  list(
    Sb = (Mb - Mw)/r,
    Sw = Mw
  )
}

#' @rdname manova
#' @param M1,M2,M3 Sum-of-squares matrices within sire, dam and individual groups.
#' @param I1,I2,I3 Number of sires, dams per sire and individuals per dam.
#' @export
manova_2way_mat <- nest_cov_estimate(manova_1way)

#' @rdname manova
#' @export
#' @param df Data frame for 2-way design.
#' @param ... Other arguments.
manova_2way_df <- df_cov_estimate(manova_2way_mat)
