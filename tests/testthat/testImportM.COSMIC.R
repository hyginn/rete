# testImportM.COSMIC.R
#
#
context("import COSMIC mutations")

# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()

testPath <- tempdir()
testFName <- paste("importM_COSMIC_", Sys.Date(), ".rds", sep = "")
# ==== END SETUP AND PREPARE ===================================================


test_that("importM.COSMIC returns error when file not found", {

    expect_error(
        importM.COSMIC(fName = "testCosmicCNA2.tsv"))
})

test_that("importM.COSMIC returns error when fNameCNA is not a valid string", {
    expect_error(
        importM.COSMIC(fNameCNA = 0))
})

test_that("importM.COSMIC returns error when silent is not a valid logical", {
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA.tsv",
                       silent = 0))
})

test_that("importM.COSMIC returns error when writeLog is not a valid logical", {
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA.tsv",
                       writeLog = 0))
})

test_that("importM.COSMIC returns error when writeToFile is not a valid logical", {
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA.tsv",
                       writeToFile = 0))
})

test_that("importM.COSMIC saves to RDS file only when writeToFile is TRUE", {
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA.tsv",
                       outFName = "cosmicCNA.rds",
                       writeToFile = FALSE))
})

test_that("importM.COSMIC only accepts outFNames of type .rds.", {
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA.tsv",
                       outFName = "cosmicCNA.txt",
                       writeToFile = TRUE))
})

test_that("importM.COSMIC returns error when outFNames is not string.", {
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA.tsv",
                       outFName = 0,
                       writeToFile = TRUE))
})

test_that("importM.COSMIC correctly returns rCNA from Cosmic CNA data", {
    fName <- "testCosmicCNA.tsv"

    expectedMatrix <- matrix(nrow = 10,
                       ncol = 1)

    expected <- data.frame(expectedMatrix)

    rownames(expected) <- c("DAZ4",
                            "FAM106A",
                            "LCE3C",
                            "PRY",
                            "RBMY1D_ENST00000382680",
                            "GYG2P1",
                            "DAZ1_ENST00000382510",
                            "ENSG00000196115",
                            "RBMY1E",
                            "DAZ2")

    colnames(expected) <- c("TCGA-683665-611825")

    expected["DAZ4", "TCGA-683665-611825"] <- 8
    expected["FAM106A", "TCGA-683665-611825"] <- -2
    expected["LCE3C", "TCGA-683665-611825"] <- -2
    expected["PRY", "TCGA-683665-611825"] <- 5
    expected["RBMY1D_ENST00000382680", "TCGA-683665-611825"] <- -2
    expected["GYG2P1", "TCGA-683665-611825"] <- 3
    expected["DAZ1_ENST00000382510", "TCGA-683665-611825"] <- 5
    expected["ENSG00000196115", "TCGA-683665-611825"] <- -2
    expected["RBMY1E", "TCGA-683665-611825"] <- -2
    expected["DAZ2", "TCGA-683665-611825"] <- 5

    expect_equal(importM.COSMIC(fName,
                                writeToFile = FALSE,
                                silent = TRUE,
                                writeLog = FALSE),
                 expected)
})

test_that("importM.COSMIC outputs rCNA RDS to file.", {
    fName <- "testCosmicCNA.tsv"
    outName <- paste(testPath, testFName, sep = "")

    expected <- importM.COSMIC(fName,
                               outName,
                               writeToFile = TRUE,
                               silent = TRUE,
                               writeLog = FALSE)

    result <- readRDS(outName)

    expect_equal(result,
                 expected)
})

test_that("importM.COSMIC returns error when CNA input data is corrupt.", {
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNAEmpty.tsv")
        )
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
