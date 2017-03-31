# testFilter
#
#
context("Tests for importFilterHypermutators.R")

# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}

SNVfileName <- 'devSNV.rds'
CNAfileName <- 'devCNA.rds'

filteredSNVfileName <- paste(getwd(), "/filteredHypermutators_", basename(SNVfileName), sep = "")
filteredCNAfileName <- paste(getwd(), "/filteredHypermutators_", basename(CNAfileName), sep = "")
# ==== END SETUP AND PREPARE ===================================================


# ==== importFilterHypermutators() =============================================

# set up a tempdir and tempfiles for outputs of importFilterHypermutators
test_that("importFilterHypermutators works correctly on valid input", {
    # run importFilterHypermutators on valid input
    expect_error(importFilterHypermutators(rSNVFileIn = c(SNVfileName),
                                           rCNAFileIn = c(CNAfileName),
                                           dOut = getwd(),
                                           xS = 400,
                                           silent = FALSE,
                                           writeLog = TRUE), NA)

    if (file.exists(filteredSNVfileName)) { file.remove(filteredSNVfileName)}
    if (file.exists(filteredCNAfileName)) { file.remove(filteredCNAfileName)}
})

test_that("importFilterHypermutators rejects invalid dOut arguments", {
    # check for NULL object
    expect_error(importFilterHypermutators(rSNVFileIn = c(SNVfileName),
                                           rCNAFileIn = c(CNAfileName),
                                           dOut = NULL,
                                           xS = 400,
                                           silent = FALSE,
                                           writeLog = TRUE))
    # check for paths that don't exist
    expect_error(importFilterHypermutators(rSNVFileIn = c(SNVfileName),
                                           rCNAFileIn = c(CNAfileName),
                                           dOut = 'no/such/path',
                                           xS = 400,
                                           silent = FALSE,
                                           writeLog = TRUE))
})

test_that("importFilterHypermutators rejects invalid xS arguments", {
    # check for including NULL objects
    expect_error(importFilterHypermutators(rSNVFileIn = c(SNVfileName),
                                           rCNAFileIn = c(CNAfileName),
                                           xS = NULL,
                                           silent = FALSE,
                                           writeLog = TRUE))
    # check character
    expect_error(importFilterHypermutators(rSNVFileIn = c(SNVfileName),
                                           rCNAFileIn = c(CNAfileName),
                                           xS = "400",
                                           silent = FALSE,
                                           writeLog = TRUE))
})

test_that("importFilterHypermutators correctly removes hypermutators", {
    # run importFilterHypermutators with a single hypermutator
    testCNA <- readRDS(CNAfileName)
    testSNV <- readRDS(SNVfileName)

    # start off with a known number of samples in each CNA and SNV file
    expect_equal(length(colnames(testCNA)[4:length(colnames(testCNA))]), 579)
    expect_equal(length(testSNV$Tumor_Sample_Barcode), 15)

    expect_error(importFilterHypermutators(rSNVFileIn = c(SNVfileName),
                                           rCNAFileIn = c(CNAfileName),
                                           xS = 15,
                                           silent = FALSE,
                                           writeLog = TRUE), NA)

    # check output file
    processedCNA <- readRDS(filteredCNAfileName)
    processedSNV <- readRDS(filteredSNVfileName)

    # expect that 458 samples are removed from rCNA (579 - 458 = 121)
    expect_equal(length(colnames(processedCNA)[4:length(colnames(processedCNA))]), 121)

    # expect that no samples are removed from rSNV
    expect_equal(length(processedSNV$Tumor_Sample_Barcode), 15)

    # test cleanup
    if (file.exists(filteredSNVfileName)) { file.remove(filteredSNVfileName)}
    if (file.exists(filteredCNAfileName)) { file.remove(filteredCNAfileName)}
})

# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)

# another check for function artifacts
if (file.exists(filteredSNVfileName)) { file.remove(filteredSNVfileName)}
if (file.exists(filteredCNAfileName)) { file.remove(filteredCNAfileName)}
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
