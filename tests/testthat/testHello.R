# testHello.R
#
#
context("Demo test")

test_that("hello.R returns a vector with the first seven Fibonacci numbers", {
    expect_equal(sum(hello()), 33)
})

# [END]