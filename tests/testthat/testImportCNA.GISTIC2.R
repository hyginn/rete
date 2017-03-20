#
# test importCNA.GISTIC2.R
#
#library(testthat)
testthat::context( "import rCNA data from GISTIC2 files")
#
#

# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()

testthat::test_that("No arguments given", {
    expect_error(
        importCNA.GISTIC2()
    )
})

testthat::test_that("Missing argument", {
    expect_error(
        importCNA.GISTIC2("inst/extdata/dCNA")
    )
})


testthat::test_that("File does not exist", {
    expect_error(
        importCNA.GISTIC2("foo.txt","inst/extdata/dCNA")
        )
})

testthat::test_that("Dir does not exist", {
    expect_error(
        importCNA.GISTIC2("tests/testthat/importhistictest.txt", "foo/")
    )
})

testthat::test_that("Not logical", {
    expect_error(
        importCNA.GISTIC2("inst/extdata/dCNA/rCNA1.rds","inst/extdata/dCNA",
                          NULL, FALSE)
    )
})

testthat::test_that("Not logical", {
    expect_error(
        importCNA.GISTIC2("inst/extdata/dCNA/rCNA1.rds","inst/extdata/dCNA",
                          FALSE, NULL)
    )
})

testthat::test_that("bad input", {
    expect_error(
        importCNA.GISTIC2("tests/testthat/importgistictest.txt","inst/extdata/dCNA",
                          FALSE, FALSE)
    )
})

# this should cover output tests
testthat::test_that("It works", {
    importCNA.GISTIC2("inst/extdata/devCNA.txt","inst/extdata/dCNA", FALSE, FALSE)
    expect_equal(readRDS("inst/extdata/dCNA/rCNA1.rds"),
        readRDS('tests/testthat/dCNA/test.rds'))
})

# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
