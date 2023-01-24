test_that("Conditional means of balanced data", {
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

  params <- list(
    list(mean = c(0, 0)),
    list(mean = c(0, 0))
  )

  balanced_data <- balance_data(y_data, params)

  expect_s3_class(balanced_data, "fullsibdata")
  expect_equal(
    dim(balanced_data$tables[[1]]),
    dim(balanced_data$tables[[2]])
  )
  expect_true(is.balanced(balanced_data))
})
