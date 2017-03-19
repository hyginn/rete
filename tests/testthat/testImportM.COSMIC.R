# testImportM.COSMIC.R
#
#
context("import COSMIC mutations")

test_that("importM.COSMIC correctly returns rCNA from Cosmic CNV data", {
    fName <- "testCosmicCNA.tsv"

    expected <- matrix(nrow = 10,
                       ncol = 1)

    rownames(expected) <- c(106410,
                            107031,
                            55218,
                            66769,
                            106400,
                            101525,
                            106281,
                            84805,
                            69781,
                            69785)

    colnames(expected) <- c("TCGA-683665-611825")

    expected["106410", "TCGA-683665-611825"] <- 10
    expected["107031", "TCGA-683665-611825"] <- 0
    expected["55218", "TCGA-683665-611825"] <- 0
    expected["66769", "TCGA-683665-611825"] <- 7
    expected["106400", "TCGA-683665-611825"] <- 0
    expected["101525", "TCGA-683665-611825"] <- 5
    expected["106281", "TCGA-683665-611825"] <- 7
    expected["84805", "TCGA-683665-611825"] <- 0
    expected["69781", "TCGA-683665-611825"] <- 0
    expected["69785", "TCGA-683665-611825"] <- 7

    expect_equal(importM.COSMIC(fName,
                                writeToFile = FALSE,
                                silent = TRUE,
                                writeLog = FALSE),
                 expected)
})

test_that("importM.COSMIC outputs rCNA RDS to file.", {
    fName <- "testCosmicCNA.tsv"
    outName <- "outCosmicCNA.rds"

    expected <- importM.COSMIC(fName,
                               outName,
                               writeToFile = TRUE,
                               silent = TRUE,
                               writeLog = FALSE)

    result <- readRDS("outCosmicCNA.rds")

    file.remove("outCosmicCNA.rds")

    expect_equal(result,
                 expected)
})

test_that("importM.COSMIC returns error when file not found", {

    expect_error(
        importM.COSMIC(fName = "testCosmicCNA2.tsv"))
})

test_that("importM.COSMIC returns error when fNameCNV is not a valid string", {
    expect_error(
        importM.COSMIC(fNameCNV = 0))
})

test_that("importM.COSMIC returns error when silent is not a valid logical", {
    expect_error(
        importM.COSMIC(fNameCNV = "testCosmicCNA.tsv",
                       silent = 0))
})

test_that("importM.COSMIC returns error when writeLog is not a valid logical", {
    expect_error(
        importM.COSMIC(fNameCNV = "testCosmicCNA.tsv",
                       writeLog = 0))
})

test_that("importM.COSMIC returns error when writeToFile is not a valid logical", {
    expect_error(
        importM.COSMIC(fNameCNV = "testCosmicCNA.tsv",
                       writeToFile = 0))
})

test_that("importM.COSMIC saves to RDS file only when writeToFile is TRUE", {
    expect_error(
        importM.COSMIC(fNameCNV = "testCosmicCNA.tsv",
                       outFName = "cosmicCNA.rds",
                       writeToFile = FALSE))
})

test_that("importM.COSMIC only accepts outFNames of type .rds.", {
    expect_error(
        importM.COSMIC(fNameCNV = "testCosmicCNA.tsv",
                       outFName = "cosmicCNA.txt",
                       writeToFile = TRUE))
})

test_that("importM.COSMIC returns error when outFNames is not string.", {
    expect_error(
        importM.COSMIC(fNameCNV = "testCosmicCNA.tsv",
                       outFName = 0,
                       writeToFile = TRUE))
})

# [END]
