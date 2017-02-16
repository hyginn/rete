# testHello.R
#
#
context("Demo test")

test_that("hello.R returns a vector with the first seven Fibonacci numbers", {
    expect_identical(hello(), c(1, 1, 2, 3, 5, 8, 13))
})

# [END]