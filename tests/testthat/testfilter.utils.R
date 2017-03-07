context("Filter utility functions")

test_that("Barcode normalized according to specs", {
    expect_equal("TCGA.A3.3307.01A.01R.0864.07",
        .filter.utils.normalizeBarcode(
        "TCGA-A3-3307-01A-01R-0864-07", separator="-"))

    expect_equal("TCGA.A3.3307.01A.01R.0864.07",
        .filter.utils.normalizeBarcode(
        "TCGA.A3.3307.01A.01R.0864.07", separator="."))

    expect_equal("TCGA.A3.3307",
        .filter.utils.normalizeBarcode(
        "tcga-a3-3307", separator='-'))

    expect_error(
        .filter.utils.normalizeBarcode(NULL))

    expect_error(
        .filter.utils.normalizeBarcode(3))

    expect_equal(
        .filter.utils.normalizeBarcode(c("TCGA-A3", "tcga-a3", "tcga-a3-3307")),
        c('TCGA.A3', 'TCGA.A3', 'TCGA.A3.3307'))
})

test_that("Patient correctly extracted from barcode", {
    expect_equal("TCGA-A3-3307",
        .filter.utils.patientFromBarcode(
        "TCGA-A3-3307-01A-01R-0864-07", separator="-"))

    expect_equal("TCGA.A3.3307",
        .filter.utils.patientFromBarcode(
        "TCGA.A3.3307.01A.01R.0864.07", separator="."))

    expect_equal("TCGA.A3.3307",
        .filter.utils.patientFromBarcode(
        "TCGA.A3.3307.01A.01R.0864.07"))

    expect_error(
        .filter.utils.patientFromBarcode("TCGA-A3"))

    expect_error(
        .filter.utils.patientFromBarcode(""))

    expect_error(
        .filter.utils.patientFromBarcode(NULL))

    expect_error(
        .filter.utils.patientFromBarcode(3))

    expect_error(
        .filter.utils.patientFromBarcode(c("TCGA", "A3", "3307")))
})

test_that("Path correctly concatenated on Linux", {
    skip_on_os(c("windows", "mac", "solaris"))

    expect_equal("/tmp/foo.rds",
        .filter.utils.outputFilename("/tmp/", "foo.rds"))

    expect_equal("/tmp/foo.rds",
        .filter.utils.outputFilename("/tmp", "foo.rds"))

    expect_equal("/tmp/foo.rds",
        .filter.utils.outputFilename("/tmp///", "foo.rds"))

    expect_equal("/tmp/foo.rds",
        .filter.utils.outputFilename("/tmp", "/home/foo/foo.rds"))
})

test_that("Path correctly concatenated on OS X", {
    skip_on_os(c("windows", "linux", "solaris"))

    expect_equal("/private/tmp/foo.rds",
        .filter.utils.outputFilename("/tmp/", "foo.rds"))

    expect_equal("/private/tmp/foo.rds",
        .filter.utils.outputFilename("/tmp", "foo.rds"))

    expect_equal("/private/tmp/foo.rds",
        .filter.utils.outputFilename("/tmp///", "foo.rds"))

    expect_equal("/private/tmp/foo.rds",
        .filter.utils.outputFilename("/tmp", "/Users/foo/foo.rds"))
})

test_that("Magic check for hash character using raw vectors", {
    maf <- writeBin(c("#version 2.2"), raw())
    gz <- writeBin(c("ZZZ"), raw())

    expect_equal("MAF",
        .filter.utils.getMagic(maf))

    expect_equal("RDS",
        .filter.utils.getMagic(gz))

    expect_error(.filter.utils.getMagic())
    expect_error(.filter.utils.getMagic(NULL))
    expect_error(.filter.utils.getMagic("not_a_real_file"))

    skip_on_os(c("windows"))
    expect_error(.filter.utils.getMagic("/dev/null"))
})

test_that("Processing an rSNV MAF file", {
    outFile <- tempfile()
    rMutName <- tempfile()
    rMutData <- c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t"),
        paste(c("TF2", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "8CD1CF38-DAD9-42EA-8B06-00955109F3D0"), collapse="\t"),
        paste(c("TF2", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0002-123-01W-0615-10",
            "E6F5767F-E87B-4FC9-8880-CA70F434EEC9"), collapse="\t"),
        paste(c("HFE", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "F813F0B2-CE88-41FB-805E-40997E0E0309"), collapse="\t"),
        paste(c("HFE", "6", "127", "456", "+", "Missense_Mutation", "SNP",
            "A", "T", "G", "TCGA-A3-0002-124-01W-0615-10", 
            "7719241D-B6C8-4B13-80F6-3047C8BBFE1F"), collapse="\t"))
    writeLines(rMutData, con=rMutName)
    expect_true(file.exists(rMutName), info="Test setup (rSNV)")
    
    filter <- c()

    output <- .filter.utils.filterrSNV(rMutName, outFile, filter)
    filtered <- readLines(outFile)
    expect_equal(filtered, c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t"),
        paste(c("TF2", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "8CD1CF38-DAD9-42EA-8B06-00955109F3D0"), collapse="\t"),
        paste(c("TF2", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0002-123-01W-0615-10",
            "E6F5767F-E87B-4FC9-8880-CA70F434EEC9"), collapse="\t"),
        paste(c("HFE", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "F813F0B2-CE88-41FB-805E-40997E0E0309"), collapse="\t"),
        paste(c("HFE", "6", "127", "456", "+", "Missense_Mutation", "SNP",
            "A", "T", "G", "TCGA-A3-0002-124-01W-0615-10", 
            "7719241D-B6C8-4B13-80F6-3047C8BBFE1F"), collapse="\t")),
        info=paste(c("All rSNV records kept with empty filter", outFile)))
    expect_equal(output[1], 4, info="Kept 4 records (empty filter)")
    expect_equal(output[2], 0, info="Removed 0 records (empty filter)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }
    
    filter <- c('ZNF660')

    output <- .filter.utils.filterrSNV(rMutName, outFile, filter)
    filtered <- readLines(outFile)
    expect_equal(filtered, c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t"),
        paste(c("TF2", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "8CD1CF38-DAD9-42EA-8B06-00955109F3D0"), collapse="\t"),
        paste(c("TF2", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0002-123-01W-0615-10",
            "E6F5767F-E87B-4FC9-8880-CA70F434EEC9"), collapse="\t"),
        paste(c("HFE", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "F813F0B2-CE88-41FB-805E-40997E0E0309"), collapse="\t"),
        paste(c("HFE", "6", "127", "456", "+", "Missense_Mutation", "SNP",
            "A", "T", "G", "TCGA-A3-0002-124-01W-0615-10", 
            "7719241D-B6C8-4B13-80F6-3047C8BBFE1F"), collapse="\t")),
        info=paste(c("All rSNV records kept", outFile)))
    expect_equal(output[1], 4, info="Kept 4 records (filter gene not in data)")
    expect_equal(output[2], 0,
        info="Removed 0 records (filter gene not in data")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }

    outFile <- tempfile()
    filter <- c('TF2')

    output <- .filter.utils.filterrSNV(rMutName, outFile, filter)
    filtered <- readLines(outFile)
    expect_equal(filtered, c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t"),
        paste(c("HFE", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-123-01W-0615-10",
            "F813F0B2-CE88-41FB-805E-40997E0E0309"), collapse="\t"),
        paste(c("HFE", "6", "127", "456", "+", "Missense_Mutation", "SNP",
            "A", "T", "G", "TCGA-A3-0002-124-01W-0615-10", 
            "7719241D-B6C8-4B13-80F6-3047C8BBFE1F"), collapse="\t")),
        info=paste(c("Some rSNV records filtered", outFile)))
    expect_equal(output[1], 2, info="Kept 2 records (filtered one gene)")
    expect_equal(output[2], 2, info="Removed 2 records (filtered one gene)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }

    outFile <- tempfile()
    filter <- c('HFE', 'TF2')

    output <- .filter.utils.filterrSNV(rMutName, outFile, filter)
    filtered <- readLines(outFile)
    expect_equal(filtered, c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t")
            ),
        info=paste(c("All rSNV records filtered out", outFile)))
    expect_equal(output[1], 0, info="Kept 0 records (filtered all genes)")
    expect_equal(output[2], 4, info="Removed 4 records (filtered all genes)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }

    outFile <- tempfile()
    rMutData <- c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t")
        )
    writeLines(rMutData, con=rMutName)
    expect_true(file.exists(rMutName), info="Test setup (rSNV)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }
    
    output <- .filter.utils.filterrSNV(rMutName, outFile, filter)
    filtered <- readLines(outFile)
    expect_equal(filtered, c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t")
            ),
        info=paste(c("Accepts bodyless rSNV files", outFile)))
    expect_equal(output[1], 0, info="Kept 0 records (none in data)")
    expect_equal(output[2], 0, info="Removed 0 records (none in data)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }

    outFile <- tempfile()
    rMutData <- c(
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t"))
    writeLines(rMutData, con=rMutName)
    expect_true(file.exists(rMutName), info="Test setup (rSNV)")
    
    expect_error(.filter.utils.filterrSNV(rMutName, outFile, filter),
        info="stop on malformed rSNV (without data)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }

    outFile <- tempfile()
    rMutData <- c(
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t"),
        paste(c("HFE", "6", "127", "456", "+", "Missense_Mutation",
            "A", "T", "G", "TCGA-A3-0001-124-01W-0615-10", 
            "7719241D-B6C8-4B13-80F6-3047C8BBFE1F"), collapse="\t"))
    writeLines(rMutData, con=rMutName)
    expect_true(file.exists(rMutName), info="Test setup (rSNV)")
    
    expect_error(.filter.utils.filterrSNV(rMutName, outFile, filter),
        info="stop on malformed rSNV (with data)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }

    outFile <- tempfile()
    rMutData <- c(
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type", "foo",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID", "extra"), collapse="\t"),
        paste(c("HFE", "6", "127", "456", "+", "Missense_Mutation", "SNP",
            "A", "T", "G", "TCGA-A3-0001-124-01W-0615-10", 'bla',
            "7719241D-B6C8-4B13-80F6-3047C8BBFE1F"), collapse="\t"))
    writeLines(rMutData, con=rMutName)
    expect_true(file.exists(rMutName), info="Test setup (rSNV)")
    
    expect_error(.filter.utils.filterrSNV(rMutName, outFile, filter),
        info="stop on malformed rSNV (incorrect record count)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }
    file.remove(rMutName)

})

test_that("Processing an rCNA RDS file", {
    outFile <- tempfile()
    rCNAName <- tempfile()
    rCNAData <- data.frame(row.names = c("HFE", "FTO", "ZNF66"),
        "TCGA-A3-0001-123" = c(0.0, 130.4, 1.5),
        "TCGA-A3-0001-124" = c(2.9, 1.3, 0.5),
        "TCGA-A3-0002-123" = c(1.0, 1.6, 3.0))
    saveRDS(rCNAData, file=rCNAName)
    expect_true(file.exists(rCNAName), info="Test setup (rCNA)")
    
    filter <- c()

    output <- .filter.utils.filterrCNA(rCNAName, outFile, filter)
    filtered <- readRDS(outFile)

    expect_equal(nrow(filtered), 3, info="3 genes present (empty filter)")
    expect_equal(ncol(filtered), 3, info="3 samples present (empty filter)")
    expect_true('HFE' %in% rownames(filtered),
        info="HFE present (empty filter)")
    expect_true('FTO' %in% rownames(filtered),
        info="FTP present (empty filter)")
    expect_true('ZNF66' %in% rownames(filtered),
        info="ZNF66 present (empty filter)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }
    
    filter <- c('TF2')

    output <- .filter.utils.filterrCNA(rCNAName, outFile, filter)
    filtered <- readRDS(outFile)

    expect_equal(nrow(filtered), 3, info="3 genes present (filter not in rCNA)")
    expect_equal(ncol(filtered), 3,
        info="3 samples present (filterd not in rCNA)")
    expect_true('HFE' %in% rownames(filtered),
        info="HFE present (filtered not in rCNA)")
    expect_true('FTO' %in% rownames(filtered),
        info="FTO present (filtered not in rCNA)")
    expect_true('ZNF66' %in% rownames(filtered),
        info="ZNF66 present (filtered not in rCNA)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }
    
    filter <- c('HFE')

    output <- .filter.utils.filterrCNA(rCNAName, outFile, filter)
    filtered <- readRDS(outFile)

    expect_equal(nrow(filtered), 2, info="2 genes present (HFE filtered)")
    expect_equal(ncol(filtered), 3, info="3 samples present (HFE filtered)")
    expect_false('HFE' %in% rownames(filtered),
        info="HFE not present (HFE filtered)")
    expect_true('FTO' %in% rownames(filtered),
        info="FTO present (HFE filtered)")
    expect_true('ZNF66' %in% rownames(filtered),
        info="ZNF66 present (HFE filtered)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }
    
    filter <- c('HFE', 'ZNF66')

    output <- .filter.utils.filterrCNA(rCNAName, outFile, filter)
    filtered <- readRDS(outFile)

    expect_equal(nrow(filtered), 1,
        info="1 gene present (HFE, ZNF66 filtered)")
    expect_equal(ncol(filtered), 3,
        info="3 samples present (HFE, ZNF66 filtered)")
    expect_false('HFE' %in% rownames(filtered),
        info="HFE not present (HFE, ZNF66 filtered)")
    expect_true('FTO' %in% rownames(filtered),
        info="FTO present (HFE, ZNF66 filtered)")
    expect_false('ZNF66' %in% rownames(filtered),
        info="ZNF66 not present (HFE, ZNF66 filtered)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }
    
    filter <- c('HFE', 'ZNF66', 'FTO')

    output <- .filter.utils.filterrCNA(rCNAName, outFile, filter)
    filtered <- readRDS(outFile)

    expect_equal(nrow(filtered), 0,
        info="No genes present (all genes filtered)")
    expect_equal(ncol(filtered), 3,
        info="3 samples present (all genes filtered")
    expect_false('HFE' %in% rownames(filtered),
        info="HFE not present (all genes filtered)")
    expect_false('FTO' %in% rownames(filtered),
        info="FTO not present (all genes filtered)")
    expect_false('ZNF66' %in% rownames(filtered),
        info="ZNF66 not present (all genes filtered)")
    if (file.exists(outFile)) {
        file.remove(outFile)
    }
    file.remove(rCNAName)
})

