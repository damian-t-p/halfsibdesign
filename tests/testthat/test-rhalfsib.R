set.seed("19940920")


q <- 2
A <- diag(rep(1, q))
mu <- rep(0, q)

I <- 2
J <- 3
K <- 4

df <- rhalfsib(mu, A, I, A, J, A, K)

test_that("Generated tibble has the right structure", {
  expect_equal(dim(df), c(q * I * J * K, 5))
  expect_setequal(names(df), c("trait", "value", "sire", "dam", "individual"))
  expect_true(is.factor(df$individual))
  expect_type(df$value, "double")
})
