test_that("Data frames are parsed properly as fullsibdata objects", {
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

  expect_s3_class(y_data, "fullsibdata")
  expect_equal(sapply(y_data$tables, nrow),
               c(2, 1))
})
