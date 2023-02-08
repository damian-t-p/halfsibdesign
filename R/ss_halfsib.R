#' @export
ss_mats <- function(data, ...) {
  UseMethod("ss_mats", data)
}

#' @export
ss_mats.halfsibdata <- function(data) {

  I <- data$dims$I
  J <- data$dims$J
  K <- data$dims$J
  
  dam_means  <- data$dam_sums / K
  sire_means <- rowsum(dam_means, data$sires) / J
  grand_mean <- colMeans(sire_means)

  ind_mom   <- data$sos
  dam_mom   <- t(dam_means) %*% dam_means
  sire_mom  <- t(sire_means) %*% sire_means
  grand_mom <- grand_mean %o% grand_mean

  list(
    m_ind  = 1 / (I * J * (K-1)) * (ind_mom - K * dam_mom),
    m_dam  = K / (I * (J - 1)) * (dam_mom - J * sire_mom),
    m_sire = (J * K) / (I - 1) * (sire_mom - I * grand_mom),
    I      = I,
    J      = J,
    K      = K,
    q      = data$dims$q
  )
}
