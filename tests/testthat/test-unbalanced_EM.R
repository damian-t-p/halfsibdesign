test_that("Paired inverse works correctly", {
  A <- matrix(c(1, 2, 2, 9), nrow = 2)
  E <- matrix(c(10, 14, 14, 20), nrow = 2)
  ns <- c(0, 1, 4, 10)

  l <- lapply(ns, \(n) solve(solve(A) + n * solve(B)))
  names(l) <- ns
  
  expect_mapequal(paired_inverse(E, A, ns), l)
})
