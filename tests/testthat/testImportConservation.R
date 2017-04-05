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
                                     outName = NULL))
    # Empty string for outName
    expect_error(.importConservation("dev.phastCons100way.wigFix",
                                     outName = ""))
})


test_that("importConservation outputs a valid file", {
    temp <- paste0(tempfile(), ".rds")
    .importConservation("dev.phastCons100way.wigFix",
                        temp, silent = TRUE)
    # Check if file exists
    #expect_equal(file.exists(testFName), TRUE)
    testOut <- readRDS(temp)
    expect_true(identical(testOut[1, ], c(1, 129)))
    expect_true(identical(testOut[6, ], c(201, 71)))
    expect_equal(length(testOut), 12)

})

test_that(".annotateConservation works on generated Chromosome files", {
    # Sample myGeneData
    myGeneData <- list()
    myGeneData[["test1"]] <- structure(list(
        sym = "TEST1",
        name = "test gene 1",
        UniProt = "A1A1A1",
        RefSeq = "NM_111111",
        assembly = "GRCh38",
        chr = "12",
        strand = "+",
        geneStart = 1,
        geneEnd = 100,
        cdsStarts = c(1L, 11L, 21L),
        cdsEnds = c(10L, 20L, 29L),
        map = function(x) {if (x>0 && x<30) { return(x) } else { return(NULL) } },
        cds = "ATG TGG GCT CAG CTC CTT CTA GGA ATG TGA",
        phastCons = NULL,
        phyloP = NULL,
        xyzMap = NULL,
        xyz = NULL),
        .Names = c("sym",
                   "name", "UniProt", "RefSeq", "assembly", "chr", "strand", "geneStart",
                   "geneEnd", "cdsStarts", "cdsEnds", "map", "cds", "phastCons",
                   "phyloP", "xyzMap", "xyz"))

    myGeneData[["test2"]] <- structure(list(
        sym = "TEST2",
        name = "test gene 2",
        UniProt = "A2A2A2",
        RefSeq = "NM_222222",
        assembly = "GRCh38",
        chr = "Y",
        strand = "-",
        geneStart = 1,
        geneEnd = 100,
        cdsStarts = c(21L, 11L, 1L),
        cdsEnds = c(29L, 20L, 10L),
        map = function(x) {if (x>0 && x<30) { return(29 - x) } else { return(NULL) } },
        cds = "ATG AAA CCC TTT GGG GGG TTT CCC AAA TGA",
        phastCons = NULL,
        phyloP = NULL,
        xyzMap = NULL,
        xyz = NULL),
        .Names = c("sym",
                   "name", "UniProt", "RefSeq", "assembly", "chr", "strand", "geneStart",
                   "geneEnd", "cdsStarts", "cdsEnds", "map", "cds", "phastCons",
                   "phyloP", "xyzMap", "xyz"))
    chromosomeScores <- readRDS("devChromosomeScores.rds")
    test <- myGeneData[["test1"]]
    scores <- .fetchConservationScores(chromosomeScores,
                                       myGeneData[["test1"]]$geneStart,
                                       myGeneData[["test1"]]$geneEnd,
                                       myGeneData[["test1"]]$cdsStarts,
                                       myGeneData[["test1"]]$cdsEnds)
    expect_equal(58, length(scores))
    expect_identical(scores[29,], c(29, 33))

})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
