# testImportConservation.R
#
#

context("Import UCSC Conservation Scores")


# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()

testPath <- tempdir()
# ==== END SETUP AND PREPARE ===================================================



test_that("parameters sent to .importConservation are correctly handled", {
    # No Parameters passed on
    expect_error(.importConservation())
    # Empty string for fName
    expect_error(.importConservation(fName = "",
                                     "testOut.rds", "phastCons"))
    # Input file is not found
    expect_error(.importConservation(fName = "dev.phastCons100way1.wigFix",
                                     "testOut.rds", "phastCons"))
    # Invalid type for fName
    expect_error(.importConservation(fName = 1, "testOut.rds", "phastCons"))
    # Invalid type for outName
    expect_error(.importConservation("dev.phastCons100way.wigFix",
                                     outName = NULL, "phastCons"))
    # Empty string for outName
    expect_error(.importConservation("dev.phastCons100way.wigFix",
                                     outName = "", "phastCons"))
    # Invalid type for scoring type
    expect_error(.importConservation("dev.phastCons100way.wigFix", "out",
                                     NULL))
    # Invalid type for silent
    expect_error(.importConservation("dev.phastCons100way.wigFix", "out",
                                     "phastCons", silent = 1))
    # Invalid type for writeLog
    expect_error(.importConservation("dev.phastCons100way.wigFix", "out",
                                     "phastCons", writeLog = NULL))
})


test_that("importConservation outputs a valid file", {
    temp <- paste0(tempfile(), ".rds")
    .importConservation("dev.phastCons100way.wigFix",
                        temp, "phastCons", silent = TRUE, writeLog = FALSE)
    # Check if file exists
    #expect_equal(file.exists(testFName), TRUE)
    testOut <- readRDS(temp)
    expect_true(identical(testOut[1, ], c(1, 129)))
    expect_true(identical(testOut[6, ], c(201, 71)))
    expect_equal(length(testOut), 12)

    # Check attached metadata
    expect_equal(attributes(testOut)$type, "ChromosomeConservationScore")
    expect_equal(attributes(testOut)$chromosome, "chrM")
    expect_equal(attributes(testOut)$scoreType, "phastCons")
    expect_equal(attributes(testOut)$version, "1.0")

    uuid_regex <- "[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
    expect_match(attributes(testOut)$UUID, uuid_regex)

})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
