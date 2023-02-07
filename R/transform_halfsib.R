#' Balance an unbalanced dataset
#'
#' @export
balance <- function(data, ...) {
  UseMethod("balance", data)
}


#' @export
balance.halfsibdata <- function(data, means, globmean = rep(0, data$dims$q)) {

  # Conditional mean of observations. These are the same within dams, so they
  # are indexed by dam
  obs_means  <- means$dam + means$sire[data$sires, ] + rep(globmean, each = nrow(means$dam))

  # Create a matrix of unobserved sums
  n_missing  <- data$dims$J - data$n.observed$inds[rownames(obs_means)]
  unobs_sums <- obs_means * n_missing

  full_sums  <- data$dam_sums + unobs_sums

  # Add unobserved sum-of-squares
  full_sos   <- data$sos + t(unobs_sums) %*% unobs_sums

  new_halfsibdata(
    sos      = full_sos,
    dam_sums = full_sums,
    sires    = data$sires,
    I        = data$dims$I,
    J        = data$dims$J,
    K        = data$dims$K,
    obs_dams = data$n.observed$dams * 0 + data$dims$J,
    obs_inds = data$n.observed$inds * 0 + data$dims$K
  )
}
