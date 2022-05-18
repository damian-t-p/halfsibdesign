set.seed(19940920)

q <- 2
A <- diag(rep(1, 2))
mu <- rep(0, 2)

I <- 2
J <- 3
K <- 4

df <- rhalfsib(mu, A, I, A, J, A, K)

expected_A <- matrix(c(16.393137, 1.330857, 1.330857, 0.108044), ncol=2)
expected_B <- matrix(c(1.858499, -1.430635, -1.430635, 1.885812), ncol=2)
expected_E <- matrix(c(0.90497485, 0.07645298, 0.07645298, 0.74779640), ncol=2)

test_that("Dimenstionality check", {
  ss <- ss_mats(df)
  expect_equal(ss$q, q)
  expect_equal(ss$I, I)
  expect_equal(ss$J, J)
  expect_equal(ss$K, K)
})

test_that("Matrix check", {
  ss <- ss_mats(df)
  expect_equal(ss$m_ind, expected_E, tolerance = 1e-5)
  expect_equal(ss$m_dam, expected_B, tolerance = 1e-5)
  expect_equal(ss$m_sire, expected_A, tolerance = 1e-5)
})
