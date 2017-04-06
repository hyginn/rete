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

# No arguments given
testthat::test_that("No arguments given", {
    expect_error(
        importCNA.GISTIC2()
    )
})

# Missing argument
testthat::test_that("Missing argument", {
    expect_error(
        importCNA.GISTIC2("inst/extdata/dCNA")
    )
})

#Source file does not exist
testthat::test_that("File does not exist", {
    expect_error(
        importCNA.GISTIC2("foo.txt","inst/extdata/dCNA")
        )
})

# Dir to write output does not exist
testthat::test_that("Dir does not exist", {
    expect_error(
        importCNA.GISTIC2("tests/testthat/importhistictest.txt", "foo/")
    )
})

# Faulty input for Boolean parameter
testthat::test_that("Not logical", {
    expect_error(
        importCNA.GISTIC2("inst/extdata/dCNA/rCNA1.rds","inst/extdata/dCNA",
                          NULL, FALSE)
    )
})


# Faulty input for Boolean parameter
testthat::test_that("Not logical", {
    expect_error(
        importCNA.GISTIC2("inst/extdata/dCNA/rCNA1.rds","inst/extdata/dCNA",
                          FALSE, NULL)
    )
})
# Invalid input file because it's only one line and all numbers for columns
testthat::test_that("bad input", {
    expect_error(
        importCNA.GISTIC2("tests/testthat/importgistictest.txt","inst/extdata/dCNA",
                          FALSE, FALSE)
    )
})

# this should cover output tests. The data in test.rds was created using
#a different method (read.delim())
testthat::test_that("It works", {
    importCNA.GISTIC2("inst/extdata/devCNA.txt","inst/extdata/dCNA", FALSE, FALSE)
    expect_equal(readRDS("inst/extdata/dCNA/rCNA1.rds"),
        readRDS('tests/testthat/dCNA/devCNA.txt.rds'))
})

# NA value in input
testthat::test_that("It works", {
    importCNA.GISTIC2("tests/testthat/devCNAwithNAvalue.txt","inst/extdata/dCNA", FALSE, FALSE)
    expect_equal(readRDS("inst/extdata/dCNA/rCNA1.rds"),
                 readRDS('tests/testthat/dCNA/rCNA1.rds'))
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
