library(testthat)

source("lib/validators.R")

test_that("returns error and quits when input value is wrong type", {
  expect_true(TRUE)

  test_that("min is not numeric", {
    validate_func <- validate_numeric(min="some string")
    expect_error(validate_func("flag", 1, list()), "'min' validate_numeric value for .+ must be numeric")
  })

  test_that("max is not numeric", {
    validate_func <- validate_numeric(max="some string")
    expect_error(validate_func("flag", 1, list()), "'max' validate_numeric value for .+ must be numeric")
  })

  test_that("min_warning error", {
    expect_true(TRUE)

    test_that("min_warning is not a list", {
      validate_func <- validate_numeric(min_warning="some string")
      expect_error(validate_func("flag", 1, list()), "'min_warning' validate_numeric value for .+ must be a length 2 list of \\(numeric, character\\)")
    })

    test_that("min_warning is wrong length", {
      validate_func <- validate_numeric(min_warning=list(1))
      expect_error(validate_func("flag", 1, list()), "'min_warning' validate_numeric value for .+ must be a length 2 list of \\(numeric, character\\)")
    })

    test_that("min_warning first value is not numeric", {
      validate_func <- validate_numeric(min_warning=list("some string", "some other string"))
      expect_error(validate_func("flag", 1, list()), "'min_warning' validate_numeric value for .+ must be a length 2 list of \\(numeric, character\\)")
    })
  })
  
  test_that("max_warning error", {
    expect_true(TRUE)

    test_that("max_warning is not a list", {
      validate_func <- validate_numeric(max_warning="some string")
      expect_error(validate_func("flag", 1, list()), "'max_warning' validate_numeric value for .+ must be a length 2 list of \\(numeric, character\\)")
    })

    test_that("max_warning is wrong length", {
      validate_func <- validate_numeric(max_warning=list(1))
      expect_error(validate_func("flag", 1, list()), "'max_warning' validate_numeric value for .+ must be a length 2 list of \\(numeric, character\\)")
    })

    test_that("max_warning first value is not numeric", {
      validate_func <- validate_numeric(max_warning=list("some string", "some other string"))
      expect_error(validate_func("flag", 1, list()), "'max_warning' validate_numeric value for .+ must be a length 2 list of \\(numeric, character\\)")
    })
  })
})

test_that("checks bounds for each flag", {
  expect_true(TRUE)
  
  test_that("min", {
    validate_func <- validate_numeric(min=1)
    expect_error(validate_func("flag", 0, list()), "Flag flag must be greater than or equal to 1 \\(current value: 0\\)")
  })

  test_that("max", {
    validate_func <- validate_numeric(max=1)
    expect_error(validate_func("flag", 2, list()), "Flag flag must be less than or equal to 1 \\(current value: 2\\)")
  })

  test_that("min_warning", {
    validate_func <- validate_numeric(min_warning=list(1, "extra info"))
    expect_warning(validate_func("flag", 0, list()), "Warning: for best results, flag flag should be greater than or equal to 1 \\(current value: 0\\); extra info")
  })

  test_that("max_warning", {
    validate_func <- validate_numeric(max_warning=list(1, "extra info"))
    expect_warning(validate_func("flag", 2, list()), "Warning: for best results, flag flag should be less than or equal to 1 \\(current value: 2\\); extra info")
  })
})

test_that("no errors or warning when all validation is met", {
  validate_func <- validate_numeric(min=0, max=4, min_warning=list(1, "too small"), max_warning=list(3, "too big"))
  expect_no_error(validate_func("flag", 2, list()))
})