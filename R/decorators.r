nest_cov_estimate <- function(f) {
  function(M1, I1, M2, I2, M3, I3) {
    S <- f(M1, I1, M2, I2)

    M_new <- S$Sw + I1 * S$Sb
    I_new <- I1 * I2

    out <- list(
      S1 = S$Sw,
      S2 = S$Sb,
      S3 = f(M_new, I_new, M3, I3)$Sb
    )

    out$reml_crit <- reml_crit(M1, out$S1, I1, M2, out$S2, I2, M3, out$S3, I3)

    return(out)
  }
}

df_cov_estimate <- function(f) {
  function(df, ...) {
    ss <- ss_mats(df)
    f(ss$m_ind, ss$K, ss$m_dam, ss$J, ss$m_sire, ss$I, ...)
  }
}
