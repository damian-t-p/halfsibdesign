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

cmean_reml <- cond_mean_reml(prior_covs, ccov_reml, data)
cmean      <- cond_mean(prior_covs, ccov, data,
                        prior_mean = cmean_reml$grand)

