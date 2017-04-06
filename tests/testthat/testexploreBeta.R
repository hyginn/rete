# testexploreBeta.R
#
#
context("explore Beta testing")

test_that("exploreBeta returns valid answer", {
    expect_true(exploreBeta(gG) <= 1 && exploreBeta(gG) >= 0)
})

# [END]