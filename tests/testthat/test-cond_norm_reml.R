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

ccov     <- cond_cov(prior_covs, data)
ccov_raw <- cond_cov_counts(prior_covs, data)

ccov_reml <- cond_cov_reml(prior_covs, ccov, ccov_raw, data)

Omega <- solve(full_prec(prior_covs, data, method = "REML"))
idx   <- cov_idx(data)

dam_names <- split(names(data$sires), data$sires)

test_that("ccov_reml is correct", {
  expect_equal(Omega[idx("group"), idx("group")], ccov_reml("group", "group"))
  for(sire in names(data$n.observed$dams)) {
    expect_equal(Omega[idx("group"), idx(sire, "group")], ccov_reml("group", sire, NA, "group"))
    expect_equal(Omega[idx(sire, "group"), idx("group")], ccov_reml(sire, "group", "group", NA))
    expect_equal(Omega[idx(sire, "group"), idx(sire, "group")], ccov_reml(sire, sire, "group", "group"))

    for(dam in dam_names[[sire]]) {
      expect_equal(Omega[idx("group"), idx(sire, dam)], ccov_reml("group", sire, NA, dam))
      expect_equal(Omega[idx(sire, dam), idx("group")], ccov_reml(sire, "group", dam, NA))
    }
    
    for(sire2 in names(data$n.observed$dams)) {
      expect_equal(Omega[idx(sire, "group"), idx(sire2, "group")],
                   ccov_reml(sire, sire2, "group", "group"))
      expect_equal(Omega[idx(sire2, "group"), idx(sire, "group")],
                   ccov_reml(sire2, sire, "group", "group"))

      for(dam in dam_names[[sire]]) {
        expect_equal(Omega[idx(sire, dam), idx(sire2, "group")],
                     ccov_reml(sire, sire2, dam, "group"))
        expect_equal(Omega[idx(sire2, "group"), idx(sire, dam)],
                     ccov_reml(sire2, sire, "group", dam))

        for(dam2 in dam_names[[sire2]]) {
          expect_equal(Omega[idx(sire, dam), idx(sire2, dam2)],
                       ccov_reml(sire, sire2, dam, dam2))
        }
      }
    }
  }
})

cmean_reml <- cond_mean_reml(prior_covs, ccov_reml, data)
cmean      <- cond_mean(prior_covs, ccov, data,
                        prior_mean = cmean_reml$grand)



