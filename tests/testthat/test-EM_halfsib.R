set.seed(12345)

q <- 4

I <- 20
J <- 3
K <- 5

mu <- 1:q

X <- rnorm(q*q) %>% matrix(nrow = q)
Q <- eigen(X, symmetric = TRUE)$vectors

A <- diag(c(3, 2, 0, 0))
B <- Q %*% diag(c(5, 0, 0, 1)) %*% t(Q)
E <- diag(q)

df_full <- rhalfsib(mu, A, I, B, J, E, K)

df_unbalanced <- df_full %>%
  dplyr::group_split(sire, dam, individual) %>%
  sample(size = floor(0.5 * (I*J*K)), replace = FALSE) %>%
  dplyr::bind_rows()

data <- halfsibdata(df_unbalanced, ind_name = individual)

prior_covs <- list(
  ind  = diag(q),
  dam  = diag(q),
  sire = diag(q)
)

EM_mats <- EM_fit.halfsibdata(data, prior_covs)

EM_mats_reml <- EM_fit.halfsibdata(
  data,
  method = "REML"
  ## err.tol = 0.01
)

###

suppressWarnings({

  lme1 <- nlme::lme(
    fixed       = value ~ -1 + trait,
    data        = df_unbalanced,
    random      = ~ -1 + trait | sire/dam,
    weights     = nlme::varIdent(form = ~ 1 | trait),
    correlation = nlme::corSymm(form =  ~ 1 | sire/dam/individual),
    method      = "REML",
    control     = nlme::lmeControl(opt = "optim", returnObject = TRUE)
  )

})

get_covs <- function(fit){
  
  resid_sds <- fit$sigma * sqrt((1 + c(0, fit$modelStruct$varStruct)))
  resid_corr <- as.matrix(fit$modelStruct$corStruct)[[1]]
  S1 <- outer(resid_sds, resid_sds) * resid_corr
  
  S2 <- fit$modelStruct$reStruct$dam %>%
      as.matrix()  %>%
      {. * fit$sigma^2}
  rownames(S2) <- NULL
  colnames(S2) <- NULL
  
  S3 <- fit$modelStruct$reStruct$sire %>%
      as.matrix()  %>%
      {. * fit$sigma^2}
  rownames(S3) <- NULL
  colnames(S3) <- NULL
  
  out <- list(
    ind  = S1,
    dam  = S2,
    sire = S3
  )
  
  return(out)
}

lme_mats <- get_covs(lme1)

test_that("EM and lme give approximately the same results", {
  expect_true(all(abs(lme_mats$ind - EM_mats$ind) < 0.1))
  expect_true(all(abs(lme_mats$dam - EM_mats$dam) < 0.1))
  expect_true(all(abs(lme_mats$sire - EM_mats$sire) < 0.1))
})


