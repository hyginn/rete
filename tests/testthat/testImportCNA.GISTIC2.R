#
# test importCNA.GISTIC2.R
#
library(testthat)
context( "import rCNA data from GISTIC2 files")
#
#


test_that("File does not exist", {
    expect_error(
        importCNA.GISTIC2(NULL,"inst/extdata/dCNA")
        )
})

test_that("Not logical", {
    expect_error(
        importCNA.GISTIC2("inst/extdata/dCNA/rCNA1.rds","inst/extdata/dCNA",
                          NULL, FALSE)
    )
})

test_that("Not logical", {
    expect_error(
        importCNA.GISTIC2("inst/extdata/dCNA/rCNA1.rds","inst/extdata/dCNA",
                          FALSE, NULL)
    )
})

test_that("bad input", {
    expect_error(
        importCNA.GISTIC2("tests/testthat/importhistictest.txt","inst/extdata/dCNA",
                          FALSE, FALSE)
    )
})


test_that("It works", {
    importCNA.GISTIC2("inst/extdata/devCNA.txt","inst/extdata/dCNA")
    expect_equal(readRDS("inst/extdata/dCNA/rCNA1.rds"),
        readRDS('tests/testthat/dCNA/testimport.rds'))
})



