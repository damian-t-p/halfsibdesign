set.seed(1234)

sires <- c(1, 1, 1, 2, 2, 3, 3, 3)
dams  <- c(1, 1, 2, 3, 4, 5, 6, 6)

q <- 3
n <- length(sires)
values <- matrix(rnorm(n * q), ncol = q)

data <- halfsibdata(values, sires, dams)

test_that("Dimensions recorded correctly", {
  expect_mapequal(
    data$dims,
    list(q = q, I = 3, J = 2, K = 2)
  )
})

test_that("Dam-sire reference works as expected", {
  expect_mapequal(
    data$sires,
    setNames(
      make.names(c(1, 1, 2, 2, 3, 3)),
      nm = make.names(c(1, 2, 3, 4, 5, 6))
    )
  )
})

test_that("Observation counts are correct", {

  expect_mapequal(
    data$n.observed$dams,
    setNames(
      c(2, 2, 2),
      nm = make.names(c(1, 2, 3))
    )
  )

  expect_mapequal(
    data$n.observed$inds,
    setNames(
      c(2, 1, 1, 1, 1, 2),
      nm = make.names(c(1, 2, 3, 4, 5, 6))
    )
  )
})
