# testMultiply.R
#
#
context("Demo test")

test_that("multiply() multiplies", {
    expect_equal(multiply(19, 3), 57)
    expect_equal(multiply(1:2, 3:4), c(3, 8))
})

# [END]