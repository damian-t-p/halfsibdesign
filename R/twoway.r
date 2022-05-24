#' Compute covariance estimates for half-sib balanced designs
#'
#' Generic function sfor computing MANOVA estimates for a balanced half-sib design
#' using various algorithms
#'
#' @param object An object for which a sum-of-squares matrix can be computed.
#' @param ... Other arguments.
#'
#' @seealso \code{\link{manova_2way_mat}}, \code{\link{preml_2way_mat}},
#' \code{\link{stepreml_2way_mat}} for computing estimates from sum-of-squares
#' matrices; \code{\link{manova_2way_df}}, \code{\link{preml_2way_df}},
#' \code{\link{stepreml_2way_df}} from a table of data.
#'
#' @name twoway

#' @rdname twoway
#' @export
manova_2way <- function(object, ...) UseMethod("manova_2way")

setMethod("manova_2way", signature=c(object="matrix"),
          function(object, ...) manova_2way_mat(object, ...)
)

setMethod("manova_2way", signature=c(object="data.frame"),
          function(object, ...) manova_2way_df(object, ...)
)

#' @rdname twoway
#' @export
preml_2way <- function(object, ...) UseMethod("preml_2way")

setMethod("preml_2way", signature=c(object="matrix"),
          function(object, ...) preml_2way_mat(object, ...)
)

setMethod("preml_2way", signature=c(object="data.frame"),
          function(object, ...) preml_2way_df(object, ...)
)


#' @rdname twoway
#' @export
stepreml_2way <- function(object, ...) UseMethod("stepreml_2way")

setMethod("stepreml_2way", signature=c(object="matrix"),
          function(object, ...) stepreml_2way_mat(object, ...)
)

setMethod("stepreml_2way", signature=c(object="data.frame"),
          function(object, ...) stepreml_2way_df(object, ...)
)
