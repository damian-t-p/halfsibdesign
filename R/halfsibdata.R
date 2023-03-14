#' Create a new object of class `halfsibdata`
#' 
#' @export
halfsibdata <- function(y, ...) {
  UseMethod("halfsibdata", y)
}

#' @export
halfsibdata.data.frame <- function(df,
                                   value_name = value,
                                   trait_name = trait,
                                   ind_name   = ind,
                                   sire_name  = sire,
                                   dam_name   = dam,
                                   df_format  = c("long", "wide"),
                                   ...) {

  df_format <- match.arg(df_format)

  if(df_format == "long") {
    df_wide <- tidyr::pivot_wider(
      df,
      id_cols     = c({{sire_name}}, {{dam_name}}, {{ind_name}}),
      names_from  = {{trait_name}},
      values_from = {{value_name}}
    )
  } else {
    df_wide <- df
  }

  value_matrix <- df_wide %>%
    dplyr::select(-c({{sire_name}}, {{dam_name}}, {{ind_name}})) %>%
    as.matrix()

  if(any(is.na(value_matrix))) {
    stop("Data cannot have missing traits")
  }
  
  halfsibdata(
    value_matrix,
    dplyr::pull(df_wide, {{sire_name}}),
    dplyr::pull(df_wide, {{dam_name}}),
    level_names  = rlang::enquos(
      sire_name  = sire_name,
      dam_name   = dam_name,
      ind_name   = ind_name,
      trait_name = trait_name,
      value_name = value_name
    ),
    ...
  )
}

#' @export
halfsibdata.matrix <- function(values, sires, dams,
                               level_names = list(
                                 sire_name  = rlang::sym("sire"),
                                 dam_name   = rlang::sym("dam"),
                                 ind_name   = rlang::sym("ind"),
                                 trait_name = rlang::sym("trait"),
                                 value_name = rlang::sym("value")
                               ),
                               ...) {
  validate_halfsibdata(values, sires, dams, level_names, ...)
}

validate_halfsibdata <- function(values, sires, dams, level_names, allow_unbalanced = TRUE) {

  sires <- make.names(sires)
  dams  <- make.names(dams)
  
  n <- nrow(values)

  if(n != length(sires) | n != length(dams)) {
    stop(
      "`sires` and `dams` must have the same lengths as the number of rows in `values`",
      call. = FALSE
    )
  }
  
  observed_dams <- tapply(
    dams,
    sires,
    \(v) length(unique(v))
  )

  if(isFALSE(allow_unbalanced) && length(unique(observed_dams)) != 1) {
    stop(
      "Must have the same number of dams for each sire",
      call. = FALSE
    )
  }
  
  observed_inds <- tapply(
    1:n,
    dams,
    \(v) length(v)
  )

  unique_dam_idx     <- !duplicated(dams)
  sire_of_dam        <- sires[unique_dam_idx]
  names(sire_of_dam) <- dams[unique_dam_idx]
  
  new_halfsibdata(
    sos         = t(values) %*% values,
    dam_sums    = rowsum(values, group = dams),
    sires       = sire_of_dam[names(observed_inds)],
    I           = length(unique(sires)),
    J           = max(observed_dams),
    K           = max(observed_inds),
    obs_dams    = c(observed_dams),
    obs_inds    = c(observed_inds),
    level_names = level_names
  )
  
}

new_halfsibdata <- function(sos,
                            dam_sums,
                            sires,
                            I, J, K,
                            obs_dams, obs_inds,
                            level_names) {

  stopifnot(is.matrix(sos))
  stopifnot(is.matrix(dam_sums))
  stopifnot(is.vector(sires))
  
  structure(
    list(
      level_names = level_names,
      sos         = sos,
      dam_sums    = dam_sums,
      sires       = sires,
      dims        = list(
        q = ncol(sos),
        I = I,
        J = J,
        K = K
      ),
      n.observed = list(
        dams = setNames(as.integer(obs_dams), names(obs_dams)),
        inds = setNames(as.integer(obs_inds), names(obs_inds))
      )
    ),
    class = "halfsibdata"
  )
  
}

is.dam_balanced <- function(data) {
  length(unique(data$n.observed$dams)) == 1
}

## is.balanced <- function(data) {
##   isTRUE(is.dam_balanced) && (length(unique(data$n.observed$inds)) == 1)
## }
