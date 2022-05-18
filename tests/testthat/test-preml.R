set.seed(19940920)

q <- 4 # number of traits

I <- 100 # number of sires
J <- 3 # number of dams
K <- 5 # number of individuals per line

mu <- 1:q

sigma_a <- 5
sigma_b <- 3
sigma_e <- 1

A <- sigma_a^2 * diag(q)
B <- sigma_b^2 * diag(q)
E <- sigma_e^2 * diag(q)

df <- rhalfsib(mu, A, I, B, J, E, K)
ss <- ss_mats(df)

test_that("Check MANOVA estimates", {
  fit_ss <- preml_2way(ss$m_ind, ss$K, ss$m_dam, ss$J, ss$m_sire, ss$I)
  fit_df <- preml_2way_df(df)

  expect_mapequal(fit_ss, fit_df)
})
