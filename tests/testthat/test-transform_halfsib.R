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
cmean <- cond_mean(prior_covs, ccov, data)

balanced_data <- balance(data, cmean)

test_that("Fully-observed observations are unchanged", {
  expect_equal(
    balanced_data$dam_sums[c("X1", "X6"), ],
    data$dam_sums[c("X1", "X6"), ]
  )
})

test_that("Unobserved observations are changed", {
  expect_false(
    any(
      balanced_data$dam_sums[c("X2", "X3"), ] ==
      data$dam_sums[c("X2", "X3"), ]
    )
  )
})

test_that("Balancing is idempotent", {
  expect_equal(
    balanced_data,
    balance(balanced_data, cmean)
  )
})

set.seed(1234)

sires <- c(1, 1, 1, 2, 2, 3, 3, 3)
dams  <- c(1, 1, 2, 3, 3, 5, 6, 6)

q <- 3
n <- length(sires)
values <- matrix(rnorm(n * q), ncol = q)

unbalanced_data <- halfsibdata(values, sires, dams)
dam_unbalanced_data <- include_unobs_dams(unbalanced_data)

ccov <- cond_cov(prior_covs, dam_unbalanced_data)
cmean <- cond_mean(prior_covs, ccov, dam_unbalanced_data)

balanced_data <- balance(dam_unbalanced_data, cmean)
