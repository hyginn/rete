context("Filter utility functions")

test_that("Patient correctly extracted from barcode", {
    expect_equal("TCGA-A3-3307",
        filter.utils.patientFromBarcode(
        "TCGA-A3-3307-01A-01R-0864-07"))

    expect_equal("TCGA.A3.3307",
        filter.utils.patientFromBarcode(
        "TCGA.A3.3307.01A.01R.0864.07", separator="."))

    expect_error(
        filter.utils.patientFromBarcode("TCGA-A3"))

    expect_error(
        filter.utils.patientFromBarcode(""))

    expect_error(
        filter.utils.patientFromBarcode(NULL))

    expect_error(
        filter.utils.patientFromBarcode(3))

    expect_error(
        filter.utils.patientFromBarcode(c("TCGA", "A3", "3307")))
})

test_that("Path correctly concatenated on Linux", {
    skip_on_os(c("windows", "mac", "solaris"))

    expect_equal("/tmp/foo.rds",
        filter.utils.outputFilename("/tmp/", "foo.rds"))

    expect_equal("/tmp/foo.rds",
        filter.utils.outputFilename("/tmp", "foo.rds"))

    expect_equal("/tmp/foo.rds",
        filter.utils.outputFilename("/tmp///", "foo.rds"))

    expect_equal("/tmp/foo.rds",
        filter.utils.outputFilename("/tmp", "/home/foo/foo.rds"))
})

test_that("Path correctly concatenated on OS X", {
    skip_on_os(c("windows", "linux", "solaris"))

    expect_equal("/private/tmp/foo.rds",
        filter.utils.outputFilename("/tmp/", "foo.rds"))

    expect_equal("/private/tmp/foo.rds",
        filter.utils.outputFilename("/tmp", "foo.rds"))

    expect_equal("/private/tmp/foo.rds",
        filter.utils.outputFilename("/tmp///", "foo.rds"))

    expect_equal("/private/tmp/foo.rds",
        filter.utils.outputFilename("/tmp", "/Users/foo/foo.rds"))
})

test_that("Magic check for hash character using raw vectors", {
    maf <- writeBin(c("#version 2.2"), raw())
    gz <- writeBin(c("ZZZ"), raw())

    expect_equal("MAF",
        filter.utils.getMagic(maf))

    expect_equal("RDS",
        filter.utils.getMagic(gz))

    expect_error(filter.utils.getMagic())
    expect_error(filter.utils.getMagic(NULL))
    expect_error(filter.utils.getMagic("not_a_real_file"))

    skip_on_os(c("windows"))
    expect_error(filter.utils.getMagic("/dev/null"))
})

test_that("Processing an rMUT MAF file", {
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
    expect_true(file.exists(rMutName), info="Test setup (rMUT)")
    
    filter <- c()

    output <- filter.utils.processrMUT(rMutName, outFile, filter)
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
        info=paste(c("All rMUT records kept with empty filter", outFile)))
    expect_equal(output[1], 4)
    expect_equal(output[2], 0)

    
    filter <- c('ZNF660')

    output <- filter.utils.processrMUT(rMutName, outFile, filter)
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
        info=paste(c("All rMUT records kept", outFile)))
    expect_equal(output[1], 4)
    expect_equal(output[2], 0)

    outFile <- tempfile()
    filter <- c('TF2')

    output <- filter.utils.processrMUT(rMutName, outFile, filter)
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
        info=paste(c("Some rMUT records filtered", outFile)))
    expect_equal(output[1], 2)
    expect_equal(output[2], 2)

    outFile <- tempfile()
    filter <- c('HFE', 'TF2')

    output <- filter.utils.processrMUT(rMutName, outFile, filter)
    filtered <- readLines(outFile)
    expect_equal(filtered, c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t")
            ),
        info=paste(c("All rMUT records filtered out", outFile)))
    expect_equal(output[1], 0)
    expect_equal(output[2], 4)

    outFile <- tempfile()
    rMutData <- c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t")
        )
    writeLines(rMutData, con=rMutName)
    expect_true(file.exists(rMutName), info="Test setup (rMUT)")
    
    output <- filter.utils.processrMUT(rMutName, outFile, filter)
    filtered <- readLines(output)
    expect_equal(filtered, c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t")
            ),
        info=paste(c("Accepts bodyless rMUT files", outFile)))
    expect_equal(output[1], 0)
    expect_equal(output[2], 0)

    outFile <- tempfile()
    rMutData <- c(
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t"))
    writeLines(rMutData, con=rMutName)
    expect_true(file.exists(rMutName), info="Test setup (rMUT)")
    
    expect_error(filter.utils.processrMUT(rMutName, outFile, filter))

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
    expect_true(file.exists(rMutName), info="Test setup (rMUT)")
    
    expect_error(filter.utils.processrMUT(rMutName, outFile, filter))

    outFile <- tempfile()
    rMutData <- c(
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID", "extra"), collapse="\t"),
        paste(c("HFE", "6", "127", "456", "+", "Missense_Mutation", "SNP",
            "A", "T", "G", "TCGA-A3-0001-124-01W-0615-10", 
            "7719241D-B6C8-4B13-80F6-3047C8BBFE1F"), collapse="\t"))
    writeLines(rMutData, con=rMutName)
    expect_true(file.exists(rMutName), info="Test setup (rMUT)")
    
    expect_error(filter.utils.processrMUT(rMutName, outFile, filter))

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

    output <- .filter.utils.processrCNA(rCNAName, outName, filter)
    filtered <- readRDS(outName)

    expect_equal(nrow(filtered), 3)
    expect_equal(ncol(filtered), 2)
    expect_false('TCGA.A3.0001.124' %in% colnames(filtered))
})

