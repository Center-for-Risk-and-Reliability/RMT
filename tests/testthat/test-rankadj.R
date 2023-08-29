test_that("rank adjustment with duplicates is correct", {
  expect_equal(rankadj(c(70,71,75,78,78,89), c(80,80,84)), c(1, 2, 3, 4, 7))
})

test_that("rank adjustment without duplicates is correct", {
  expect_equal(rankadj(c(70,71,75,78,89), c(80,80,84)), c(1, 2, 3, 4, 6.5))
})