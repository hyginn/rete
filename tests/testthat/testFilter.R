# testFilter
#
#
context("Tests for importFilterHypermutators.R")

# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
# ==== END SETUP AND PREPARE ===================================================


# ==== importFilterHypermutators() =============================================

# set up a tempdir and tempfiles for outputs of importFilterHypermutators
test_that("importFilterHypermutators works correctly on valid input", {
    # run importFilterHypermutators on valid input
    expect_error(importFilterHypermutators(rSNVFileIn = c('../../inst/extdata/devSNV.rds'),
                                           rCNAFileIn = c('../../inst/extdata/devCNA.rds'),
                                           dOut = getwd(),
                                           xS = 400,
                                           silent = FALSE,
                                           writeLog = TRUE), NA)
})

test_that("importFilterHypermutators rejects invalid dOut arguments", {
    # check for NULL object
    expect_error(importFilterHypermutators(rSNVFileIn = c('../../inst/extdata/devSNV.rds'),
                                           rCNAFileIn = c('../../inst/extdata/devCNA.rds'),
                                           dOut = NULL,
                                           xS = 400,
                                           silent = FALSE,
                                           writeLog = TRUE))
    # check for paths that don't exist
    expect_error(importFilterHypermutators(rSNVFileIn = c('../../inst/extdata/devSNV.rds'),
                                           rCNAFileIn = c('../../inst/extdata/devCNA.rds'),
                                           dOut = 'no/such/path',
                                           xS = 400,
                                           silent = FALSE,
                                           writeLog = TRUE))
})

test_that("importFilterHypermutators rejects invalid xS arguments", {
    # check for including NULL objects
    expect_error(importFilterHypermutators(rSNVFileIn = c('../../inst/extdata/devSNV.rds'),
                                           rCNAFileIn = c('../../inst/extdata/devCNA.rds'),
                                           xS = NULL,
                                           silent = FALSE,
                                           writeLog = TRUE))
    # check character
    expect_error(importFilterHypermutators(rSNVFileIn = c('../../inst/extdata/devSNV.rds'),
                                           rCNAFileIn = c('../../inst/extdata/devCNA.rds'),
                                           xS = "400",
                                           silent = FALSE,
                                           writeLog = TRUE))
})

test_that("importFilterHypermutators correctly removes hypermutators", {
    # run importFilterHypermutators with a single hypermutator

    # check output file
})

# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
