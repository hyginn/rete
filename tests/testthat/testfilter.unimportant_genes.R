context("Unimportant gene filtering")

test_that("Patient correctly extracted from barcode", {
    expect_equal("TCGA-A3-3307",
        .filter.unimportant_genes.patientFromBarcode(
        "TCGA-A3-3307-01A-01R-0864-07"))

    expect_equal("TCGA.A3.3307",
        .filter.unimportant_genes.patientFromBarcode(
        "TCGA.A3.3307.01A.01R.0864.07", separator="."))

    expect_error(
        .filter.unimportant_genes.patientFromBarcode("TCGA-A3"))

    expect_error(
        .filter.unimportant_genes.patientFromBarcode(""))

    expect_error(
        .filter.unimportant_genes.patientFromBarcode(NULL))

    expect_error(
        .filter.unimportant_genes.patientFromBarcode(3))

    expect_error(
        .filter.unimportant_genes.patientFromBarcode(c("TCGA", "A3", "3307")))
})

test_that("Barcode correctly extracted from hashkey", {
    expect_equal("TCGA-A3-3307",
        .filter.unimportant_genes.barcodeFromHashkey(
        "HFE|TCGA-A3-3307"))

    expect_equal("TCGA.A3.3307",
        .filter.unimportant_genes.barcodeFromHashkey(
        "ZNF660^TCGA.A3.3307", separator="^"))

    expect_error(
        .filter.unimportant_genes.barcodeFromHashkey("|TCGA-A3"))

    expect_error(
        .filter.unimportant_genes.barcodeFromHashkey("TCGA-A3"))

    expect_error(
        .filter.unimportant_genes.barcodeFromHashkey(""))

    expect_error(
        .filter.unimportant_genes.barcodeFromHashkey(NULL))

    expect_error(
        .filter.unimportant_genes.barcodeFromHashkey(3))

    expect_error(
        .filter.unimportant_genes.barcodeFromHashkey(c("TCGA", "A3", "3307")))
})

test_that("Path correctly concatenated on Linux", {
    skip_on_os(c("windows", "mac", "solaris"))

    expect_equal("/tmp/foo.rds",
        .filter.unimportant_genes.outputFilename("/tmp/", "foo.rds"))

    expect_equal("/tmp/foo.rds",
        .filter.unimportant_genes.outputFilename("/tmp", "foo.rds"))

    expect_equal("/tmp/foo.rds",
        .filter.unimportant_genes.outputFilename("/tmp///", "foo.rds"))

    expect_equal("/tmp/foo.rds",
        .filter.unimportant_genes.outputFilename("/tmp", "/home/foo/foo.rds"))
})

test_that("Path correctly concatenated on OS X", {
    skip_on_os(c("windows", "linux", "solaris"))

    expect_equal("/private/tmp/foo.rds",
        .filter.unimportant_genes.outputFilename("/tmp/", "foo.rds"))

    expect_equal("/private/tmp/foo.rds",
        .filter.unimportant_genes.outputFilename("/tmp", "foo.rds"))

    expect_equal("/private/tmp/foo.rds",
        .filter.unimportant_genes.outputFilename("/tmp///", "foo.rds"))

    expect_equal("/private/tmp/foo.rds",
        .filter.unimportant_genes.outputFilename("/tmp", "/Users/foo/foo.rds"))
})

test_that("Magic check for hash character using raw vectors", {
    maf <- writeBin(c("#version 2.2"), raw())
    gz <- writeBin(c("ZZZ"), raw())

    expect_equal("MAF",
        .filter.unimportant_genes.getMagic(maf))

    expect_equal("RDS",
        .filter.unimportant_genes.getMagic(gz))

    expect_error(.filter.unimportant_genes.getMagic())
    expect_error(.filter.unimportant_genes.getMagic(NULL))
    expect_error(.filter.unimportant_genes.getMagic("not_a_real_file"))

    skip_on_os(c("windows"))
    expect_error(.filter.unimportant_genes.getMagic("/dev/null"))
})

test_that("Loading expression rCNA data", {
    rCNApath <- tempfile()
    listy <- data.frame(row.names = c("HFE", "FTO", "ZNF66"),
        a = c(0.0, 130.4, 1.5),
        b = c(2.9, 1.3, 0.5),
        c = c(1.0, 1.6, 0.0))
    expect_is(listy, "data.frame", info="Test setup (raw data)")

    saveRDS(listy, file=rCNApath)
    expect_true(file.exists(rCNApath), info="Test setup (RDS)")

    expression <- .filter.unimportant_genes.loadExpression(rCNApath, 2)
    expect_equal(expression$numSampleReads, 2)

    expect_equal(expression$samples$"HFE|a", list(ignore=FALSE, numReads=0.0))
    expect_equal(expression$samples$"FTO|a", list(ignore=FALSE, numReads=130.4))
    expect_equal(expression$samples$"ZNF66|a", list(ignore=FALSE, numReads=1.5))

    expect_equal(expression$samples$"HFE|b", list(ignore=FALSE, numReads=2.9))
    expect_equal(expression$samples$"FTO|b", list(ignore=FALSE, numReads=1.3))
    expect_equal(expression$samples$"ZNF66|b", list(ignore=FALSE, numReads=0.5))

    expect_equal(expression$samples$"HFE|c", list(ignore=FALSE, numReads=1.0))
    expect_equal(expression$samples$"FTO|c", list(ignore=FALSE, numReads=1.6))
    expect_equal(expression$samples$"ZNF66|c", list(ignore=FALSE, numReads=0.0))

    expression <- .filter.unimportant_genes.loadExpression(rCNApath, 500)
    expect_equal(expression$numSampleReads, 0)

    expect_equal(expression$samples$"HFE|a", list(ignore=FALSE, numReads=0.0))
    expect_equal(expression$samples$"FTO|a", list(ignore=FALSE, numReads=130.4))
    expect_equal(expression$samples$"ZNF66|a", list(ignore=FALSE, numReads=1.5))

    expect_equal(expression$samples$"HFE|b", list(ignore=FALSE, numReads=2.9))
    expect_equal(expression$samples$"FTO|b", list(ignore=FALSE, numReads=1.3))
    expect_equal(expression$samples$"ZNF66|b", list(ignore=FALSE, numReads=0.5))

    expect_equal(expression$samples$"HFE|c", list(ignore=FALSE, numReads=1.0))
    expect_equal(expression$samples$"FTO|c", list(ignore=FALSE, numReads=1.6))
    expect_equal(expression$samples$"ZNF66|c", list(ignore=FALSE, numReads=0.0))

    expression <- .filter.unimportant_genes.loadExpression(rCNApath, 1)
    expect_equal(expression$numSampleReads, 5)

    expect_equal(expression$samples$"HFE|a", list(ignore=FALSE, numReads=0.0))
    expect_equal(expression$samples$"FTO|a", list(ignore=FALSE, numReads=130.4))
    expect_equal(expression$samples$"ZNF66|a", list(ignore=FALSE, numReads=1.5))

    expect_equal(expression$samples$"HFE|b", list(ignore=FALSE, numReads=2.9))
    expect_equal(expression$samples$"FTO|b", list(ignore=FALSE, numReads=1.3))
    expect_equal(expression$samples$"ZNF66|b", list(ignore=FALSE, numReads=0.5))

    expect_equal(expression$samples$"HFE|c", list(ignore=FALSE, numReads=1.0))
    expect_equal(expression$samples$"FTO|c", list(ignore=FALSE, numReads=1.6))
    expect_equal(expression$samples$"ZNF66|c", list(ignore=FALSE, numReads=0.0))

    writeLines(text=c("blah", "blah"), con=rCNApath)
    expect_error(.filter.unimportant_genes.loadExpression(rCNApath, 1))
})

test_that("Filtering expression data", {
    rCNApath <- tempfile()
    listy <- data.frame(row.names = c("HFE", "FTO", "ZNF66"),
        "TCGA-A3-0001-123" = c(0.0, 130.4, 1.5),
        "TCGA-A3-0001-124" = c(2.9, 1.3, 0.5),
        "TCGA-A3-0002-123" = c(1.0, 1.6, 3.0))
    expect_is(listy, "data.frame", info="Test setup (raw data)")

    saveRDS(listy, file=rCNApath)
    expect_true(file.exists(rCNApath), info="Test setup (RDS)")

    expression <- .filter.unimportant_genes.loadExpression(rCNApath, 2)
    expect_equal(expression$numSampleReads, 2)
   
    vital <- new.env(hash=TRUE)
    vital[["TCGA.A3.0001"]] <- 5
    vital[["TCGA.A3.0002"]] <- 500

    filter <- .filter.unimportant_genes.filterExpression(
        expression, vital, rT=2, pT=0.75, cT=10)

    expect_equal(filter$"HFE|TCGA.A3.0001.123",
        list(ignore=FALSE, numReads=0))
    expect_equal(filter$"HFE|TCGA.A3.0001.124",
        list(ignore=FALSE, numReads=2.9))
    expect_equal(filter$"HFE|TCGA.A3.0002.123",
        list(ignore=TRUE, numReads=1.0))

    expect_equal(filter$"FTO|TCGA.A3.0001.123",
        list(ignore=FALSE, numReads=130.4))
    expect_equal(filter$"FTO|TCGA.A3.0001.124",
        list(ignore=FALSE, numReads=1.3))
    expect_equal(filter$"FTO|TCGA.A3.0002.123",
        list(ignore=TRUE, numReads=1.6))

    expect_equal(filter$"ZNF66|TCGA.A3.0001.123",
        list(ignore=FALSE, numReads=1.5))
    expect_equal(filter$"ZNF66|TCGA.A3.0001.124",
        list(ignore=FALSE, numReads=0.5))
    expect_equal(filter$"ZNF66|TCGA.A3.0002.123",
        list(ignore=FALSE, numReads=3.0))

})
