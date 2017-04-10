#
# test importCNA.GISTIC2.R
#
#library(testthat)
########## ToDo (bs):  remove comment
########## ToDo (bs):  don't use testthat::
testthat::context( "import rCNA data from GISTIC2 files")
#
#

########## ToDo (bs):  _inst_ directory does not exist in installed package
########## ToDo (bs):  change all writes to tempdir()

########## ToDo (bs):  Move test-data reads to inst subdirectory and properly
########## ToDo (bs):  fetch system file






# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()

# No arguments given
testthat::test_that("No arguments given", {
    expect_error(
        importCNA.GISTIC2()
    )
########## ToDo (bs):  Check error message
})

# Missing argument
testthat::test_that("Missing argument", {
########## ToDo (bs):  Missing _what_ argument?
    expect_error(
        importCNA.GISTIC2("inst/extdata/dCNA")
########## ToDo (bs):  Check error message
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
########## ToDo (bs):  Check error message. This should actually generate the wrong error 
########## ToDo (bs):  since you wrote a type into the input file name
    )
})

# Faulty input for Boolean parameter
testthat::test_that("Not logical", {
########## ToDo (bs):  "not logical" _what_?
    expect_error(
        importCNA.GISTIC2("inst/extdata/dCNA/rCNA1.rds","inst/extdata/dCNA",
                          NULL, FALSE)
########## ToDo (bs):  Misleading: this does not specify a valid _input_ file;
########## ToDo (bs):  ... might fail for the wrong reason.
########## ToDo (bs):  Better not to _assume_ argument order in a test.
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
########## ToDo (bs):  silent should be TRUE in testing
    )
})

# this should cover output tests. The data in test.rds was created using
#a different method (read.delim())
testthat::test_that("It works", {
########## ToDo (bs):  "It works" is not a good descriptor for a test
    importCNA.GISTIC2("inst/extdata/devCNA.txt","inst/extdata/dCNA", FALSE, FALSE)
########## ToDo (bs):  This file was not among the committed assets
########## ToDo (bs):  Don't write output there!
    expect_equal(readRDS("inst/extdata/dCNA/rCNA1.rds"),
        readRDS('tests/testthat/dCNA/devCNA.txt.rds'))
########## ToDo (bs):  Is this correct usage of the return value of readRDS()?
    
})

# NA value in input
testthat::test_that("It works", {
    importCNA.GISTIC2("tests/testthat/devCNAwithNAvalue.txt","inst/extdata/dCNA", FALSE, FALSE)
########## ToDo (bs):  CRITICAL ERROR: you are not testing the output that this 
########## ToDo (bs):  function call has produced, but the one you created before.

    expect_equal(readRDS("inst/extdata/dCNA/rCNA1.rds"),
                 readRDS('tests/testthat/dCNA/rCNA1.rds'))
})

########## ToDo (bs):  Missing: log file tests
########## ToDo (bs):  Missing: check metadata
########## ToDo (bs):  Missing: confirm function of silent and writeLog parameters
########## ToDo (bs):  Missing: test for incomplete input file
########## ToDo (bs):  How shoudld NA values _really_ get handled?



# ==== BEGIN TEARDOWN AND RESTORE ==============================================
########## ToDo (bs):  CRITICAL ERROR - test output was not removed.

logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
