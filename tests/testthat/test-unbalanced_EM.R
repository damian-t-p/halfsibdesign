test_that("Paired inverse computation", {
  A <- matrix(c(1, 2, 2, 9), nrow = 2)
  E <- matrix(c(10, 14, 14, 20), nrow = 2)
  ns <- c(0, 1, 4, 10)

  l <- lapply(ns, \(n) solve(solve(A) + n * solve(E)))
  names(l) <- ns
  
  expect_mapequal(paired_inverse(E, A, ns), l)
})

set.seed(123)

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

test_that("EM algorithm", {

  ests <- EM_oneway(y_data, diag(2), diag(2))
  
  expect_true(
    all(sapply(
      ests,
      \(X) all(abs(X - t(X)) < 1e-6) & all(eigen(X)$values >= 0)
    ))
  )

})
