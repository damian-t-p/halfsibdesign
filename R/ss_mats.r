#' Compute sum-of-squares matrices from half-sib design data
#'
#' Given a dataframe of half-sib responses \eqn{y_{ijk}}, computes the following
#' matrices:
#' \deqn{S_A = \frac{JK}{I-1}\sum_{i} (\bar{y}_{i} - \bar{y}) (y_{i} - \bar{y})'}
#' \deqn{S_B = \frac{K}{I(J-1)}\sum_{i,j} (\bar{y}_{ij} - \bar{y}_{i})
#' (\bar{y}_{ij} - \bar{y}_{i})'}
#' \deqn{S_E = \frac{1}{IJ(K-1)}\sum_{i,j,k} (y_{ijk} - \bar{y}_{ij})
#' (y_{ijk} - \bar{y}_{ij})'}
#' These matrices satisfy
#' \deqn{E S_E = E,
#' E S_B = E + K B,
#' E S_A = E + K B + JK A,}
#' where A, B, E are the population covariance matrices of the sires, dams and
#' individuals respectively.
#'
#' @param df Half-sib data frame with column names `trait`, `value`, `sire`,
#' `dam`, `individual`
#'
#' @return A list with the following entries:
#' \itemize{
#'   \item `m_ind`, `m_dam`, `m_sire`: The matrices S_E, S_B and S_A respectively
#'   \item `q`: The dimensionality of the responses
#'   \item `I`, `J`, `K`: The number of sires, dams per sire and individuals per dam
#' }
#' @export
ss_mats <- function(df) {

  q <- df$trait %>% unique %>% length
  I <- df$sire %>% unique %>% length
  J <- (df$dam %>% unique %>% length) / I
  K <- (df$individual %>% unique %>% length) / (I*J)

  df_means <- df %>%
    dplyr::group_by(trait) %>%
    dplyr::mutate(grand = mean(value)) %>%
    dplyr::group_by(trait, sire) %>%
    dplyr::mutate(sire = mean(value)) %>%
    dplyr::group_by(trait, sire, dam) %>%
    dplyr::mutate(dam = mean(value)) %>%
    dplyr::ungroup()

  vals <- df_means$value %>% matrix(nrow = q)
  grand_mean <- df_means$grand %>% matrix(nrow = q)
  sire_mean <- df_means$sire %>% matrix(nrow = q)
  dam_mean <- df_means$dam %>% matrix(nrow = q)

  list(
    m_ind = (vals - dam_mean) %*% t(vals - dam_mean) / (I*J*(K-1)),
    m_dam = (dam_mean - sire_mean) %*% t(dam_mean - sire_mean) / (I*(J-1)),
    m_sire = (sire_mean - grand_mean) %*% t(sire_mean - grand_mean) / (I-1),
    I = I,
    J = J,
    K = K,
    q = q
  )
}
