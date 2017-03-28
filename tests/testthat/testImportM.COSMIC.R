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

    # Test to make sure an error is raised when the file specified at fNameCNA
    # is not found.
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA2.tsv"))
})

test_that("importM.COSMIC returns error when fNameCNA is not a valid string", {

    # Test to make sure an error is raised when the file name specified at
    # fNameCNA is not a valid string.
    expect_error(
        importM.COSMIC(fNameCNA = 0))
})

test_that("importM.COSMIC returns error when silent is not a valid logical", {

    # Test to make sure an error is raised when silent is not a valid logical.
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA.tsv",
                       silent = 0))
})

test_that("importM.COSMIC returns error when writeLog is not a valid logical", {

    # Test to make sure an error is raised when writeLog is not a valid logical.
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA.tsv",
                       writeLog = 0))
})

test_that("importM.COSMIC returns error when outFName is not string.", {

    # Test to make sure an error is raised when outFName is not a valid string.
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA.tsv",
                       outFName = 0))
})

# Test to make sure that importM.COSMIC outputs to file specified at outFName.
test_that("importM.COSMIC outputs to file specified at outFName.", {
    fNameCNA <- "testCosmicCNA.tsv"
    outName <- paste(testPath, testFName, sep = "")

    importM.COSMIC(fNameCNA,
                   outName,
                   silent = TRUE,
                   writeLog = FALSE)

    expect_equal(file.exists(outName),
                 TRUE)
})

# Test to make sure that the data frame outputted to file is valid.
test_that("importM.COSMIC outputs correct rCNA data frame to file.", {

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

    fNameCNA <- "testCosmicCNA.tsv"
    outName <- paste(testPath, testFName, sep = "")

    importM.COSMIC(fNameCNA,
                   outName,
                   silent = TRUE,
                   writeLog = FALSE)

    result <- readRDS(outName)

    meta <- list(type = attributes(result)$type,
                 version = attributes(result)$version,
                 UUID = attributes(result)$UUID)

    for (name in names(meta)) {
        attr(expected, name) <- meta[[name]]
    }

    expect_equal(result,
                 expected)
})

# Test that attributes of rCNA data frame are correct.
test_that("importM.COSMIC assigns appropriate attributes to rCNA data frame.", {
    fNameCNA <- "testCosmicCNA.tsv"
    outName <- paste(testPath, testFName, sep = "")

    importM.COSMIC(fNameCNA,
                   outName,
                   silent = TRUE,
                   writeLog = FALSE)

    result <- readRDS(outName)

    expect_equal(attributes(result)$type, "rCNA")
    expect_equal(attributes(result)$version, "1.0")

    uuid_regex <- "[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
    expect_match(attributes(result)$UUID, uuid_regex)

})

# Tests to make sure corrupted CNA input raises an error.
test_that("importM.COSMIC returns error when CNA input data is corrupt.", {

    # Test to make sure error is raised when CNA input has negative TOTAL_CN
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA_totalCN_negative.tsv")
    )

    # Test to make sure error is raised when CNA input file does not contain
    # correct columns
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA_columns.tsv")
    )

    # Test to make sure error is raised when CNA input file has no CNA data
    # (only headings)
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNAEmpty.tsv")
        )

    # Test to make sure error is raised when CNA input file has gene_name
    # equal to "NA"
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA_gene_name_has_NA.tsv")
    )

    # Test to make sure error is raised when CNA input file has gene_name
    # equal to "N/A"
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA_gene_name_has_N_A.tsv")
    )

    # Test to make sure error is raised when CNA input file has gene_name
    # equal to "-"
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA_gene_name_has_dash.tsv")
    )

    # Test to make sure error is raised when CNA input file has ID_SAMPLE
    # of type string
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA_sample_string.tsv")
    )

    # Test to make sure error is raised when CNA input file has ID_TUMOUR
    # of type string
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA_tumour_string.tsv")
    )

    # Test to make sure error is raised when CNA input file has TOTAL_CN
    # of type string
    expect_error(
        importM.COSMIC(fNameCNA = "testCosmicCNA_totalCN_string.tsv")
    )

})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
