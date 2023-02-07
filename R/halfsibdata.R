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
                                   dam_name   = dam) {
  
  df_wide <- tidyr::pivot_wider(
    df,
    id_cols     = c({{sire_name}}, {{dam_name}}, {{ind_name}}),
    names_from  = {{trait_name}},
    values_from = {{value_name}}
  )

  value_matrix <- df_wide %>%
    dplyr::select(-c({{sire_name}}, {{dam_name}}, {{ind_name}})) %>%
    as.matrix()

  if(any(is.na(value_matrix))) {
    stop("Data cannot have missing traits")
  }
  
  halfsibdata(
    value_matrix,
    dplyr::pull(df_wide, {{sire_name}}),
    dplyr::pull(df_wide, {{dam_name}})
  )
}

#' @export
halfsibdata.matrix <- function(values, sires, dams) {
  validate_halfsibdata(values, sires, dams)
}

validate_halfsibdata <- function(values, sires, dams) {

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

  if(length(unique(observed_dams)) != 1) {
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
    sos      = t(values) %*% values,
    dam_sums = rowsum(values, group = dams),
    sires    = sire_of_dam,
    I        = length(unique(sires)),
    J        = max(observed_dams),
    K        = max(observed_inds),
    obs_dams = c(observed_dams),
    obs_inds = c(observed_inds)
  )
  
}

new_halfsibdata <- function(sos,
                            dam_sums,
                            sires,
                            I, J, K,
                            obs_dams, obs_inds) {

  stopifnot(is.matrix(sos))
  stopifnot(is.matrix(dam_sums))
  stopifnot(is.vector(sires))
  
  structure(
    list(
      sos      = sos,
      dam_sums = dam_sums,
      sires    = sires,
      dims     = list(
        q = ncol(sos),
        I = I,
        J = J,
        K = K
      ),
      n.observed = list(
        dams = obs_dams,
        inds = obs_inds
      )
    ),
    class = "halfsibdata"
  )
  
}
