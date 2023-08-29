test_that("rank calculation without censored data works", {
  expect_equal(rankcalc( c(50, 70, 88, 112, 140, 160)), c(1,2,3,4,5,6))
})

test_that("rank calculation with censored data works", {
  expect_equal(rankcalc( c(245, 312, 409),c(500, 500, 500)), c(1,2,3))
})