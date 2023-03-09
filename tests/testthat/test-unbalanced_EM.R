test_that("Paired inverse computation", {
  A <- matrix(c(1, 2, 2, 9), nrow = 2)
  E <- matrix(c(10, 14, 14, 20), nrow = 2)
  ns <- c(0, 1, 4, 10)

  l <- lapply(ns, \(n) solve(solve(A) + n * solve(E)))
  names(l) <- ns
  
  expect_mapequal(paired_inverse(E, A, ns), l)
})

set.seed(1234)

df <- tibble::tribble(
  ~sire, ~ind, ~trait,
  1,     1,    1,
  1,     1,    2,
  1,     2,    1,
  1,     2,    2,
  2,     1,    2,
  2,     1,    1
)

df$value <- rnorm(nrow(df))

y_data <- fullsibdata(
  df,
  sire_name  = sire,
  ind_name   = ind,
  trait_name = trait,
  value_name = value
)

test_that("Conditional sire effect", {

  l <- alpha_cond_params(y_data, diag(2), diag(2))
  expect_type(l, "list")
  expect_length(l, length(unique(df$sire)))
  expect_equal(names(l[[1]]), c("mean", "cov"))

})

test_that("REML EM algorithm", {

  ests <- EM_oneway(y_data, diag(2), diag(2), method = "REML")
  
  expect_true(
    all(sapply(
      ests,
      \(X) all(abs(X - t(X)) < 1e-6)
    ))
  )

  expect_true(
    all(sapply(
      ests,
      \(X) all(eigen(X)$values > -1e-6)
    ))
  )
})

test_that("ML EM algorithm", {
  ests <- EM_oneway(y_data, diag(2), diag(2), method = "ML")
  
  expect_true(
    all(sapply(
      ests,
      \(X) all(abs(X - t(X)) < 1e-6)
    ))
  )

  expect_true(
    all(sapply(
      ests,
      \(X) all(eigen(X)$values > -1e-6)
    ))
  )

})


set.seed(12345)

q <- 4

J <- 20
K <- 5

mu <- 1:q

X <- rnorm(q*q) %>% matrix(nrow = q)
Q <- eigen(X, symmetric = TRUE)$vectors

B <- Q %*% diag(c(5, 0, 0, 1)) %*% t(Q)
E <- diag(q)

df_full <- rfullsib(mu, B, J, E, K)

df_unbalanced <- df_full %>%
  dplyr::group_split(sire, individual) %>%
  sample(size = floor(0.75 * (J*K)), replace = FALSE) %>%
  dplyr::bind_rows()


data <- df_unbalanced %>%
  mutate(intercept = 1) %>%
  halfsibdata(
    sire_name = intercept,
    dam_name  = sire,
    ind_name  = individual,
  )

y_data <- fullsibdata(
  df_unbalanced,
  sire_name  = sire,
  ind_name   = individual,
  trait_name = trait,
  value_name = value
)

REML_mats <- EM_fit(data, method = "ML_nofix")
ML_mats <- EM_oneway(y_data, diag(q), diag(q), method = "ML")

test_that("REML and MLE give approximately the same results", {
  expect_true(all(abs(REML_mats$individual - ML_mats$Sigma_E) < 0.2))
  expect_true(all(abs(REML_mats$sire - ML_mats$Sigma_A) < 0.2))
})
