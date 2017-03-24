# testCheckGeneSymbols.R
#
#

context("checkGeneSymbols functions")

test_that("isGeneSymbol identifies invalid input", {
    # Input contains invalid input.
    noSuchObject <- c(TRUE, 3.14, data.frame(key=c(1,2)), c(1,2,3))
    expect_error(isGeneSymbol())
    expect_error(isGeneSymbol(c()))
    expect_error(isGeneSymbol(NULL))
    expect_error(isGeneSymbol(noSuchObject))
    expect_error(isGeneSymbol(TRUE))
    expect_error(isGeneSymbol(3.14))
    expect_error(isGeneSymbol(data.frame(key=c(1,2))))
    expect_error(isGeneSymbol(c(1,2,3)))
    expect_error(isGeneSymbol(NA))
})

test_that("isGeneSymbol identifies correct gene symbols,
          and input is case insensitive", {
    # Input contains correct HGNC gene symbols.
    expect_equal(isGeneSymbol(c("A1BG", "a1bg", "A1bG", "a1Bg")),
                              c(TRUE, TRUE, TRUE, TRUE))
})

test_that("isGeneSymbol identifies incorrect gene symbols
          that are not present in the table", {
    # Input includes symbols not present in the table.
    expect_equal(isGeneSymbol(c("1234", "!@#$@", NA, "A1BG")),
                              c(FALSE, FALSE, FALSE, TRUE))
})

# [END]
