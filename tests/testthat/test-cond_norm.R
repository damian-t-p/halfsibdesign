set.seed(1234)

sires <- c(1, 1, 1, 2, 2, 3, 3, 3)
dams  <- c(1, 1, 2, 3, 4, 5, 6, 6)

q <- 3
n <- length(sires)
values <- matrix(rnorm(n * q), ncol = q)

data <- halfsibdata(values, sires, dams)

prior_covs <- list(
  ind  = rnorm(q*q) %>% matrix(nrow = 3) %>% {. %*% t(.)},
  dam  = rnorm(q*q) %>% matrix(nrow = 3) %>% {. %*% t(.)},
  sire = rnorm(q*q) %>% matrix(nrow = 3) %>% {. %*% t(.)}
)

ccov <- cond_cov(prior_covs, data)

test_that("Covariance matrices are symmetric", {
  expect_equal(ccov("X1", "X1", "X2"), ccov("X1", "X2", "X1"))
  expect_equal(ccov("X2", "X3", "X4"), ccov("X2", "X3", "X4"))
  expect_equal(ccov("X3", "X5", "X6"), ccov("X3", "X6", "X5"))
})

test_that("Covariance matrices work correctly for dams with equal observation numbers", {
  expect_equal(data$n.observed$inds[["X3"]], data$n.observed$inds[["X4"]])
  expect_equal(ccov("X2", "X3", "X3"), ccov("X2", "X4", "X4"))
  expect_false(all(ccov("X2", "X3", "X3") == ccov("X2", "X3", "X4")))
})

cmean <- cond_mean(ccov, data)

test_that("Sire means have correct dimensions", {
  expect_equal(ncol(cmean$sire), q)
  expect_equal(rownames(cmean$sire), names(data$n.observed$dams))
})

test_that("Dam means have correct dimensions", {
  expect_equal(ncol(cmean$dams), q)
  expect_equal(rownames(cmean$dam), names(data$n.observed$inds))
})
