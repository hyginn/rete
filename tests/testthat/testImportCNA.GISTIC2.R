#
# test importCNA.GISTIC2.R
#
context( "import rCNA data from GISTIC2 files")
#
#
v <- c(NULL)
v2 <- c("inst/extdata/devCNA.txt")

test_that("File does not exist", {
    expect_error(
        importCNA.GISTIC2(v,"inst/extdata/dCNA")
        )
})

test_that("It works", {
    importCNA.GISTIC2(v2,"inst/extdata/dCNA")
    expect_equal(readRDS("inst/extdata/dCNA/rCNA1.rds"),
        readRDS('tests/testthat/dCNA/testimport.rds'))
})

