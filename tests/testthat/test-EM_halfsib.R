set.seed(1234)

q <- 4

I <- 20
J <- 3
K <- 5

mu <- 1:q

A <- diag(c(1, 1, 0, 0))
B <- diag(c(1, 0, 0, 1))
E <- diag(q)

df_full <- rhalfsib(mu, A, I, B, J, E, K)

df_unbalanced <- df_full %>%
  dplyr::group_split(sire, dam, individual) %>%
  sample(size = floor(0.9 * (I*J*K)), replace = FALSE) %>%
  dplyr::bind_rows()

data <- halfsibdata(df_unbalanced, ind_name = individual)

prior_covs <- list(
  ind  = diag(q),
  dam  = diag(q),
  sire = diag(q)
)

EM_mats <- EM_fit.halfsibdata(data, prior_covs)

###

lme1 <- nlme::lme(
  fixed       = value ~ -1 + trait,
  data        = df_unbalanced,
  random      = ~ -1 + trait | sire/dam,
  weights     = nlme::varIdent(form = ~ 1 | trait),
  correlation = nlme::corSymm(form =  ~ 1 | sire/dam/individual),
  method      = "REML",
  control     = nlme::lmeControl(opt = "optim", returnObject = TRUE)
)

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
    S1 = S1,
    S2 = S2,
    S3 = S3
  )
  
  return(out)
}

lme_mats <- get_covs(lme1)

test_that("EM and lme give approximately the same results", {
  expect_true(all(abs(lme_mats$S1 - EM_mats$ind) < 0.1))
  expect_true(all(abs(lme_mats$S2 - EM_mats$dam) < 0.1))
  expect_true(all(abs(lme_mats$S3 - EM_mats$sire) < 0.1))
})
