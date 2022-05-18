#' Compute 1- or 2-way MANOVA estimates
#'
#' Given either a tuple of sum-of-squares matrices or a dataframe, computes
#' MANOVA estimates for the covariances
#' @export
manova_1way <- function(Mw, r, Mb, n) {
  list(
    Sb = (Mb - Mw)/r,
    Sw = Mw
  )
}

#' @export
manova_2way <- nest_cov_estimate(manova_1way)

#' @export
manova_2way_df <- df_cov_estimate(manova_2way)
