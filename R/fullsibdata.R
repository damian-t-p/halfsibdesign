#' Create a new object of class `fullsibdata`
#'
#' @export
fullsibdata <- function(y, ...) {
  UseMethod("fullsibdata", y)
}

#' @export
fullsibdata.data.frame <- function(y_df,
                                   sire_name = sire,
                                   ind_name = ind,
                                   trait_name = trait,
                                   value_name = value) {
  
  sire  <- dplyr::enquo(sire_name)
  ind   <- dplyr::enquo(ind_name)
  trait <- dplyr::enquo(trait_name)
  value <- dplyr::enquo(value_name)

  split_dfs <- y_df %>%
    dplyr::arrange(!!trait) %>%
    tidyr::pivot_wider(names_from = !!trait,
                       values_from = !!value) %>%
    dplyr::group_split(!!sire)

  lapply(split_dfs, \(df) as.matrix(df[ , 3:ncol(df)])) %>%
    fullsibdata.list()
}

#' @export
fullsibdata.list <- function(y_tables) {
  validate_fullsibdata(
    new_fullsibdata(y_tables)
  )
}

validate_fullsibdata <- function(y_data) {
  stopifnot(is.list(y_data$tables))
  stopifnot(all(sapply(y_data$tables, is.matrix)))
  
  y_data
}

new_fullsibdata <- function(y_tables) {

  structure(
    list(
      tables = y_tables
    ),
    class = "fullsibdata"
  )
}

#' @export
is.propensitymodel <- function(x) inherits(x, "propensitymodel")
