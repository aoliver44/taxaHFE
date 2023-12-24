library(testthat)
source("tree.R")

test_that("returns true if greater than 3", {
    expect_true(greater_than_3(4))
})


test_that("returns false if less than 3", {
    expect_false(greater_than_3(2))
})
