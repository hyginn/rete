# testDummy.R

context("<anuvar tests>")


# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()
# ==== END SETUP AND PREPARE ===================================================



test_that("parameter errors are correctly handled", {
    # no input
    expect_error(anuvar())
    # numeric input
    expect_error(anuvar(numeric(), testgeneData))
    expect_error(anuvar(testSNV, numeric()))

    # missing geneData table
    expect_error(anuvar(SNV = testSNV))
    # missing SNV table
    expect_error(anuvar(geneData = testgeneData))
    # null input
    expect_error(anuvar(NULL, testgeneData))
    expect_error(anuvar(testSNV, NULL))

    # vector of more than size 1 input
    expect_error(anuvar(c(testSNV1, testSNV2), testgeneData))

    # invalid input for silent and writeLog
    expect_error(anuvar(SNV = testSNV, geneData = testgeneData, silent = NULL))
    expect_error(anuvar(SNV = testSNV, geneData = testgeneData, writeLog = NULL))
    expect_error(anuvar(SNV = testSNV, geneData = testgeneData, silent = 123))
    expect_error(anuvar(SNV = testSNV, geneData = testgeneData, writeLog = 123))
    expect_error(anuvar(SNV = testSNV, geneData = testgeneData, silent = c(TRUE, TRUE)))
    expect_error(anuvar(SNV = testSNV, geneData = testgeneData, writeLog = c(TRUE, TRUE)))
})


test_that("a sane input gives an expected output", {
    # make small input SNV table
    # make small input geneData table of annotations
    # make small output TRUE table for that input SNV and geneData tables
    # run anuvar module on input SNV and geneData
    # compare TRUE ouput with module OUTput
})


test_that("a corrupt input does not lead to corrupted output", {
    # check if missing required column in SNV table does not give corrupted output
    # check if missing required column form geneData able does not give corrupted output

    # check if missing values in geneData table are not added to SNV
})



# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
