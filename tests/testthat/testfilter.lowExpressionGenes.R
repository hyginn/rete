context("Low expression gene filtering")

test_that("Parsing valid expression data", {
    exprName <- tempfile()
    exprData <- c(
        paste(c("Hybridization REF",
            "TCGA-A0-0001-01A-11R-ABCD-01",
            "TCGA-A0-0001-01B-11R-ABCD-01",
            "TCGA-A0-0002-01A-11R-1234-99"), collapse="\t"),
        paste(c("gene_id", "normalized_count",
            "normalized_count", "normalized_count"), collapse="\t"),
        paste(c("HFE|3077", "173.9", "12.5", "1.0"), collapse="\t"),
        paste(c("TF2|123", "4003.1", "0.0", "6.5"), collapse="\t")
        )
    writeLines(exprData, con=exprName)
    expect_true(file.exists(exprName), info="Test setup (expression data)")

    data <- .filter.lowExpressionGenes.parseExpressionFile(exprName)
    expect_equal(nrow(data), 2)
    expect_equal(ncol(data), 3)
    expect_equal(data['TF2', 'TCGA.A0.0001.01B.11R.ABCD.01'], 0.0)
    unlink(exprName)
})

test_that("Parsing valid expression data, ignoring normals", {
    exprName <- tempfile()
    exprData <- c(
        paste(c("Hybridization REF",
            "TCGA-A0-0001-01A-11R-ABCD-01",
            "TCGA-A0-0001-01B-11R-ABCD-01",
            "TCGA-A0-0002-11A-11R-1234-99"), collapse="\t"),
        paste(c("gene_id", "normalized_count",
            "normalized_count", "normalized_count"), collapse="\t"),
        paste(c("HFE|3077", "12.5", "1.0", "7.5"), collapse="\t"),
        paste(c("TF2|123", "0.0", "6.5", "0.5"), collapse="\t")
        )
    writeLines(exprData, con=exprName)
    expect_true(file.exists(exprName), info="Test setup (expression data)")

    data <- .filter.lowExpressionGenes.parseExpressionFile(exprName)
    expect_equal(nrow(data), 2)
    expect_equal(ncol(data), 2)
    expect_equal(data['TF2', 'TCGA.A0.0001.01B.11R.ABCD.01'], 6.5)
    unlink(exprName)
})

test_that("Parsing invalid expression data", {
    exprName <- tempfile()
    exprData <- c(
        paste(c("Hybridization REF",
            "TCGA-A0-0001-01A-11R-ABCD-01",
            "TCGA-A0-0001-01B-11R-ABCD-01",
            "TCGA-A0-0002-01A-11R-1234-99"), collapse="\t"),
        paste(c("gene_id", "normalized_count",
            "normalized_count", "normalized_count"), collapse="\t"),
        paste(c("HFE|3077", "173.9", "12.5", "1.0"), collapse="\t"),
        paste(c("TF2|123", "4003.1", "foo", "6.5"), collapse="\t")
        )
    writeLines(exprData, con=exprName)
    expect_true(file.exists(exprName), info="Test setup (expression data)")

    expect_error(.filter.lowExpressionGenes.parseExpressionFile(exprName))
    unlink(exprName)
})

test_that("Parsing non-normalized RNAseq expression data", {
    exprName <- tempfile()
    exprData <- c(
        paste(c("Hybridization REF",
            "TCGA-A0-0001-01A-11R-ABCD-01",
            "TCGA-A0-0001-01A-11R-ABCD-01",
            "TCGA-A0-0001-01A-11R-ABCD-01"), collapse="\t"),
        paste(c("gene_id", 'raw_count',
            "scaled_estimate", "transcript_id"), collapse="\t"),
        paste(c("HFE|3077", "173.9", "12.5", "uc010"), collapse="\t"),
        paste(c("TF2|123", "4003.1", "0.0", "l1d"), collapse="\t")
        )
    writeLines(exprData, con=exprName)
    expect_true(file.exists(exprName), info="Test setup (expression data)")

    expect_error(.filter.lowExpressionGenes.parseExpressionFile(exprName),
        info='throws an error when run against non-normalized data')
    unlink(exprName)
})

test_that("Parsing expression data with incomplete header", {
    exprName <- tempfile()
    exprData <- c(
        paste(c("Hybridization REF",
            "TCGA-A0-0001-01A-11R-ABCD-01",
            "TCGA-A0-0001-01B-11R-ABCD-01",
            "TCGA-A0-0001-01C-11R-ABCD-01"), collapse="\t"),
        paste(c("HFE|3077", "173.9", "12.5", "22.9"), collapse="\t"),
        paste(c("TF2|123", "4003.1", "0.0", "1.0"), collapse="\t")
        )
    writeLines(exprData, con=exprName)
    expect_true(file.exists(exprName), info="Test setup (expression data)")

    expect_error(.filter.lowExpressionGenes.parseExpressionFile(exprName),
        info='throws an error when run against data that lacks both header lines')
    unlink(exprName)
})

test_that("Parsing expression data without samples", {
    exprName <- tempfile()
    exprData <- c("Hybridization REF", 'gene_id', "HFE|3077", "TF2|123")
    writeLines(exprData, con=exprName)
    expect_true(file.exists(exprName), info="Test setup (expression data)")

    expect_error(.filter.lowExpressionGenes.parseExpressionFile(exprName),
        info='throws an error when run against data that lacks samples')
    unlink(exprName)
})

test_that("Parsing empty expression data file", {
    exprName <- tempfile()
    file.create(exprName)
    expect_true(file.exists(exprName), info="Test setup (expression data)")

    expect_error(.filter.lowExpressionGenes.parseExpressionFile(exprName),
        info='throws an error when run against an empty data file')
    unlink(exprName)
})

test_that("Processing valid expression data, single file", {
    exprName <- tempfile()
    exprData <- c(
        paste(c("Hybridization REF",
            "TCGA-A0-0001-01A-11R-ABCD-01",
            "TCGA-A0-0001-01B-11R-ABCD-01",
            "TCGA-A0-0002-01A-11R-1234-99"), collapse="\t"),
        paste(c("gene_id", "normalized_count",
            "normalized_count", "normalized_count"), collapse="\t"),
        paste(c("HFE|3077", "1.0", "1.0", "1.0"), collapse="\t"),
        paste(c("TF2|123", "4003.1", "123.0", "621.5"), collapse="\t")
        )
    writeLines(exprData, con=exprName)
    expect_true(file.exists(exprName), info="Test setup (expression data)")

    data <- .filter.lowExpressionGenes.processExpressionData(exprName,
        rT=3, pT=0.7)
    expect_equal(length(data), 1)
    expect_equal(data[1], 'HFE')
    unlink(exprName)
})

test_that("Processing valid expression data, two files", {
    exprName1 <- tempfile()
    exprName2 <- tempfile()
    exprData <- c(
        paste(c("Hybridization REF",
            "TCGA-A0-0001-01A-11R-ABCD-01",
            "TCGA-A0-0001-01B-11R-ABCD-01",
            "TCGA-A0-0002-01A-11R-1234-99"), collapse="\t"),
        paste(c("gene_id", "normalized_count",
            "normalized_count", "normalized_count"), collapse="\t"),
        paste(c("HFE|3077", "1.0", "1.0", "1.0"), collapse="\t"),
        paste(c("HFE2|3077", "1.0", "1.0", "1.0"), collapse="\t"),
        paste(c("ZNF660|3077", "1.0", "1.0", "1.0"), collapse="\t"),
        paste(c("TF2|123", "4003.1", "123.0", "621.5"), collapse="\t")
        )
    writeLines(exprData, con=exprName1)
    expect_true(file.exists(exprName1), info="Test setup (expression data 1)")

    exprData <- c(
        paste(c("Hybridization REF",
            "TCGA-A0-0004-01A-11R-ABCD-01",
            "TCGA-A0-0004-01B-11R-ABCD-01",
            "TCGA-A0-0005-01A-11R-1234-99"), collapse="\t"),
        paste(c("gene_id", "normalized_count",
            "normalized_count", "normalized_count"), collapse="\t"),
        paste(c("HFE|3077", "1.0", "1.0", "1.0"), collapse="\t"),
        paste(c("HFE3|3077", "1.0", "1.0", "1.0"), collapse="\t"),
        paste(c("ZNF660|3077", "100.0", "100.0", "100.0"), collapse="\t"),
        paste(c("TF3|123", "4003.1", "123.0", "621.5"), collapse="\t")
        )
    writeLines(exprData, con=exprName2)
    expect_true(file.exists(exprName2), info="Test setup (expression data 2)")

    data <- .filter.lowExpressionGenes.processExpressionData(
        c(exprName1, exprName2),
        rT=3, pT=0.7)
    expect_equal(length(data), 4)
    expect_true('HFE' %in% data)
    expect_true('HFE2' %in% data)
    expect_true('HFE3' %in% data)
    expect_true('ZNF660' %in% data)
    unlink(exprName1)
    unlink(exprName2)
})

test_that("A full run of the main function", {
    dOut <- tempfile()
    dir.create(dOut)
    expect_true(dir.exists(dOut), info="Test setup (dOut)")

    wlName <- tempfile()
    writeLines(c('ZNF66'), con=wlName)

    exprName <- tempfile()
    exprData <- c(
        paste(c("Hybridization REF",
            "TCGA-A0-0001-01A-11R-ABCD-01",
            "TCGA-A0-0001-01B-11R-ABCD-01",
            "TCGA-A0-0002-01A-11R-1234-99"), collapse="\t"),
        paste(c("gene_id", "normalized_count",
            "normalized_count", "normalized_count"), collapse="\t"),
        paste(c("HFE|3077", "173.9", "12.5", "52.0"), collapse="\t"),
        paste(c("ZNF66|3077", "2.0", "1.5", "0.0"), collapse="\t"),
        paste(c("TF2|123", "1.1", "0.0", "0.5"), collapse="\t")
        )
    writeLines(exprData, con=exprName)
    expect_true(file.exists(exprName), info="Test setup (expression data)")

    rCNAName <- tempfile()
    rCNAData <- data.frame(row.names = c("HFE", "TF2", "ZNF66"),
        "TCGA.A3.0001.123" = c(0.0, 130.4, 1.5),
        "TCGA.A3.0001.124" = c(2.9, 1.3, 0.5),
        "TCGA.A3.0002.123" = c(1.0, 1.6, 3.0))
    saveRDS(rCNAData, file=rCNAName)
    expect_true(file.exists(rCNAName), info="Test setup (rCNA)")
    
    rSNVName <- tempfile()
    rSNVData <- c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t"),
        paste(c("ZNF66", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "F813F0B2-CE88-41FB-805E-40997E0E0309"), collapse="\t"),
        paste(c("TF2", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "F813F0B2-CE88-41FB-805E-40997E0E0309"), collapse="\t"),
        paste(c("HFE", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "F813F0B2-CE88-41FB-805E-40997E0E0309"), collapse="\t"),
        paste(c("HFE", "6", "127", "456", "+", "Missense_Mutation", "SNP",
            "A", "T", "G", "TCGA-A3-0001-124-01W-0615-10", 
            "7719241D-B6C8-4B13-80F6-3047C8BBFE1F"), collapse="\t"))
    writeLines(rSNVData, con=rSNVName)
    expect_true(file.exists(rSNVName), info="Test setup (rSNV)")

    outSNV <- .filter.utils.outputFilename(dOut, basename(rSNVName))
    outCNA <- .filter.utils.outputFilename(dOut, basename(rCNAName))

    removed <- filter.lowExpressionGenes(exprName, wlName,
        rCNA=c(rCNAName), rSNV=c(rSNVName),
        dOut, rT=3, pT=0.7, silent=TRUE, writeLog=FALSE)

    expect_equal(length(removed), 1)
    expect_equal(removed, c('TF2'))

    CNA <- readRDS(outCNA)
    expect_equal(nrow(CNA), 2)
    expect_equal(ncol(CNA), 3)
    expect_false('TF2' %in% rownames(CNA))
    expect_true('HFE' %in% rownames(CNA))
    expect_true('ZNF66' %in% rownames(CNA))
    
    SNV <- readLines(outSNV)
    keepSNV <- c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t"),
        paste(c("ZNF66", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "F813F0B2-CE88-41FB-805E-40997E0E0309"), collapse="\t"),
        paste(c("HFE", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "F813F0B2-CE88-41FB-805E-40997E0E0309"), collapse="\t"),
        paste(c("HFE", "6", "127", "456", "+", "Missense_Mutation", "SNP",
            "A", "T", "G", "TCGA-A3-0001-124-01W-0615-10", 
            "7719241D-B6C8-4B13-80F6-3047C8BBFE1F"), collapse="\t"))
    expect_equal(SNV, keepSNV)

    unlink(rCNAName)
    unlink(rSNVName)
    unlink(wlName)
    unlink(exprData)
    unlink(outSNV)
    unlink(outCNA)
    unlink(dOut)
})
