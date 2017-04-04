# testImportConservation.R
#
#

context("Import UCSC Conservation Scores")


# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()

testPath <- tempdir()
testFName <- paste("importConservation_", Sys.Date(), ".rds", sep = "")
# ==== END SETUP AND PREPARE ===================================================



test_that("parameters sent to .importConservation are correctly handled", {
    # No Parameters passed on
    expect_error(.importConservation())
    # Empty string for fName
    expect_error(.importConservation(fName = "",
                                     "testOut.rds", "M", "phastCons"))
    # Input file is not found
    expect_error(.importConservation(fName = "dev.phastCons100way1.wigFix",
                                     "testOut.rds", "M", "phastCons"))
    # Invalid type for fName
    expect_error(.importConservation(fName = 1, "testOut.rds", "M", "phastCons"))
    # Invalid type for outName
    expect_error(.importConservation("dev.phastCons100way.wigFix",
                                     outName = NULL, "M", "phastCons"))
    # Empty string for outName
    expect_error(.importConservation("dev.phastCons100way.wigFix",
                                     outName = "", "M", "phastCons"))
    # Invalid type for chromosome
    expect_error(.importConservation("dev.phastCons100way.wigFix", "out",
                                     12, "phastCons"))
    # Invalid type for scoring type
    expect_error(.importConservation("dev.phastCons100way.wigFix", "out",
                                     "M", NULL))
    # Invalid type for silent
    expect_error(.importConservation("dev.phastCons100way.wigFix", "out",
                                     "M", "phastCons", silent = 1))
    # Invalid type for writeLog
    expect_error(.importConservation("dev.phastCons100way.wigFix", "out",
                                     "M", "phastCons", writeLog = NULL))
})


test_that("importConservation outputs a valid file", {
    .importConservation("dev.phastCons100way.wigFix",
                        testFName, "1", "phastCons")
    # Check if file exists
    expect_equal(file.exists(testFName), TRUE)
    testOut <- readRDS(testFName)
    expect_true(identical(testOut[1, ], c(1, 129)))
    expect_true(identical(testOut[202, ], c(202, 71)))
    expect_equal(length(testOut), 404)

    # Check attached metadata
    expect_equal(attributes(result)$type, "Conservation Score")
    expect_equal(attributes(result)$chromosome, "1")
    expect_equal(attributes(result)$scoring, "phastCons")
    expect_equal(attributes(result)$version, "1.0")

    uuid_regex <- "[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
    expect_match(attributes(result)$UUID, uuid_regex)

})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
if (file.exists(testFName)) { file.remove(testFName)}
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
