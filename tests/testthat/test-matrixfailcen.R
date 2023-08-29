test_that("censored matrix is formed properly", {
  expect_equal(matrix.failcen(xi = c(70,71,75,78,89), rc = c(80,80,84), nx = 8), matrix(c(70, 71, 75, 78, 80, 80, 84, 89, 1, 1, 1, 1, 0, 0, 0, 1), nrow = 8, ncol = 2))
})