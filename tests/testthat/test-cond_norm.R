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
  expect_equal(ccov("X1", "X1", "X2"), t(ccov("X1", "X2", "X1")))
  expect_equal(ccov("X2", "X3", "X4"), t(ccov("X2", "X4", "X3")))
  expect_equal(ccov("X3", "X5", "X6"), t(ccov("X3", "X6", "X5")))
})

test_that("Covariance matrices work correctly for dams with equal observation numbers", {
  expect_equal(data$n.observed$inds[["X3"]], data$n.observed$inds[["X4"]])
  expect_equal(ccov("X2", "X3", "X3"), ccov("X2", "X4", "X4"))
  expect_false(all(ccov("X2", "X3", "X3") == ccov("X2", "X3", "X4")))
})

Omega <- solve(full_prec(prior_covs, data, method = "ML"))
idx   <- cov_idx(data, method = "ML")

dam_names <- split(names(data$sires), data$sires)

test_that("ccov agrees with manual inversion", {
  for(sire in names(data$n.observed$dams)) {
    expect_equal(Omega[idx(sire, "group"), idx(sire, "group")], ccov(sire, "group", "group"))
    for(dam1 in dam_names[[sire]]) {
      expect_equal(Omega[idx(sire, dam1), idx(sire, "group")], ccov(sire, dam1, "group"))
      expect_equal(Omega[idx(sire, "group"), idx(sire, dam1)], ccov(sire, "group", dam1))
      for(dam2 in dam_names[[sire]]) {
         expect_equal(Omega[idx(sire, dam1), idx(sire, dam2)], ccov(sire, dam1, dam2))
      }
    }
  }
})

cmean <- cond_mean(prior_covs, ccov, data)

test_that("Sire means have correct dimensions", {
  expect_equal(ncol(cmean$sire), q)
  expect_equal(rownames(cmean$sire), names(data$n.observed$dams))
})

test_that("Dam means have correct dimensions", {
  expect_equal(ncol(cmean$dam), q)
  expect_equal(rownames(cmean$dam), names(data$n.observed$inds))
})
