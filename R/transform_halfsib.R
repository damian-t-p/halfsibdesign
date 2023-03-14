#' Balance an unbalanced dataset
#'
#' @export
balance <- function(data, ...) {
  UseMethod("balance", data)
}


#' @export
balance.halfsibdata <- function(data, means, globmean = rep(0, data$dims$q)) {

  stopifnot(is.dam_balanced(data))
  
  # Conditional mean of observations. These are the same within dams, so they
  # are indexed by dam
  obs_means  <- means$dam + means$sire[data$sires, ] + rep(globmean, each = nrow(means$dam))

  # Create a matrix of unobserved sums
  n_missing  <- data$dims$K - data$n.observed$inds[rownames(obs_means)]
  unobs_sums <- obs_means * n_missing

  full_sums  <- data$dam_sums + unobs_sums

  # Add unobserved sum-of-squares
  unobs_sums_semi <- obs_means * sqrt(n_missing)
  full_sos   <- data$sos + t(unobs_sums_semi) %*% unobs_sums_semi

  new_halfsibdata(
    sos         = full_sos,
    dam_sums    = full_sums,
    sires       = data$sires,
    I           = data$dims$I,
    J           = data$dims$J,
    K           = data$dims$K,
    obs_dams    = data$n.observed$dams * 0L + data$dims$J,
    obs_inds    = data$n.observed$inds * 0L + data$dims$K,
    level_names = data$level_names
  )
}


include_unobs_dams <- function(data) {
  new_sires <- extend_names(
    data$dims$J - data$n.observed$dams,
    existing_dams = names(data$sires)
  )

  n_new <- length(new_sires)
  
  new_dam_sums <- matrix(
    data     = 0,
    nrow     = n_new,
    ncol     = data$dims$q,
    dimnames = list(names(new_sires))
  )

  data$dam_sums          <- rbind(data$dam_sums, new_dam_sums)
  data$sires             <- c(data$sires, new_sires)
  data$n.observed$dams[] <- data$dims$J
  data$n.observed$inds   <- c(data$n.observed$inds,
                              setNames(rep(0L, n_new), names(new_sires)))

  data
}

#' Create new names for dams based on their sires
extend_names <- function(sires, existing_dams = character(0)) {
  unaccounted_dams <- sires %>%
    mapply(\(count, sire) rep(sire, count), ., names(.)) %>%
    Reduce(c, .)

  k <- length(existing_dams)
  full_dams <- c(existing_dams, unaccounted_dams)
  
  names(unaccounted_dams) <- make.names(full_dams, unique = TRUE)[(k+1):length(full_dams)]

  unaccounted_dams
}
