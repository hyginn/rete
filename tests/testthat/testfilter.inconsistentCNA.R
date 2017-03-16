context("Filter CNA with inconsistent signs")

test_that("Filtering inconsistent CNA out of an rCNA RDS file, iT set", {
    dOut <- tempfile()
    rCNAName <- tempfile()
    filter <- tempfile()
    rCNAData <- data.frame(row.names = c("HFE", "FTO", "ZNF66", "TF2"),
        "TCGA-A3-0001-02A" = c(3.0, -130.4, -1.5, 0.0),
        "TCGA-A3-0001-02B" = c(2.9, -1.3, -0.5, 0.0),
        "TCGA-A3-0002-02A" = c(1.0, -1.6, 3.0, 0.0),
        "TCGA-A3-0002-02B" = c(1.0, -1.6, 3.0, 0.0))
    saveRDS(rCNAData, file=rCNAName)
    expect_true(file.exists(rCNAName), info="Test setup (rCNA)")
    dir.create(dOut)
    expect_true(dir.exists(dOut), info="Test setup (output directory)")
    outFile <- .filter.utils.outputFilename(dOut, rCNAName)
    
    output <- filter.inconsistentCNA(rCNAName, c(rCNAName), dOut, filter,
        iT=0.5)

    expect_equal(length(output), 1, info="1 gene filtered")
    expect_true('TF2' %in% output, info='TF2 has inconsistent sign')

    CNA <- readRDS(outFile)
    expect_equal(nrow(CNA), 3, info="3 genes present")
    expect_equal(ncol(CNA), 4, info="4 samples present")
    expect_true('HFE' %in% rownames(CNA), info="HFE present")
    expect_true('FTO' %in% rownames(CNA), info="FTO present")
    expect_true('ZNF66' %in% rownames(CNA), info="ZNF66 present")
    expect_false('TF2' %in% rownames(CNA), info="TF2 filtered out")

    genes <- readLines(con=filter)
    expect_equal(genes, output, info="output list same as filter file contents")
    
    file.remove(filter)
    file.remove(outFile)
    file.remove(rCNAName)
    unlink(dOut)
})

test_that("Filtering inconsistent CNA out of an rCNA RDS file", {
    dOut <- tempfile()
    rCNAName <- tempfile()
    filter <- tempfile()
    rCNAData <- data.frame(row.names = c("HFE", "FTO", "ZNF66", "TF2"),
        "TCGA-A3-0001-02A" = c(3.0, -130.4, -1.5, 0.0),
        "TCGA-A3-0001-02B" = c(2.9, -1.3, -0.5, 0.0),
        "TCGA-A3-0002-02A" = c(1.0, -1.6, 3.0, 0.0),
        "TCGA-A3-0002-02B" = c(1.0, -1.6, 3.0, 0.0))
    saveRDS(rCNAData, file=rCNAName)
    expect_true(file.exists(rCNAName), info="Test setup (rCNA)")
    dir.create(dOut)
    expect_true(dir.exists(dOut), info="Test setup (output directory)")
    outFile <- .filter.utils.outputFilename(dOut, rCNAName)
    
    output <- filter.inconsistentCNA(rCNAName, c(rCNAName), dOut, filter)

    expect_equal(length(output), 2, info="2 genes filtered")
    expect_true('ZNF66' %in% output, info='ZNF66 has inconsistent sign')
    expect_true('TF2' %in% output, info='TF2 has inconsistent sign')

    CNA <- readRDS(outFile)
    expect_equal(nrow(CNA), 2, info="2 genes present")
    expect_equal(ncol(CNA), 4, info="4 samples present")
    expect_true('HFE' %in% rownames(CNA), info="HFE present")
    expect_true('FTO' %in% rownames(CNA), info="FTO present")
    expect_false('ZNF66' %in% rownames(CNA), info="ZNF66 filtered out")
    expect_false('TF2' %in% rownames(CNA), info="TF2 filtered out")

    genes <- readLines(con=filter)
    expect_equal(genes, output, info="output list same as filter file contents")
    
    file.remove(filter)
    file.remove(outFile)
    file.remove(rCNAName)
    unlink(dOut)
})

test_that("Filtering inconsistent CNA out of an rSNV MAF file", {
    dOut <- tempfile()
    rCNAName <- tempfile()
    filter <- tempfile()
    rCNAData <- data.frame(row.names = c("HFE", "FTO", "ZNF66", "TF2"),
        "TCGA-A3-0001-02A" = c(3.0, -130.4, -1.5, 0.0),
        "TCGA-A3-0001-02B" = c(2.9, -1.3, -0.5, 0.0),
        "TCGA-A3-0002-02A" = c(1.0, -1.6, 3.0, 0.0),
        "TCGA-A3-0002-02B" = c(1.0, -1.6, 3.0, 0.0))
    saveRDS(rCNAData, file=rCNAName)
    expect_true(file.exists(rCNAName), info="Test setup (rCNA)")
    dir.create(dOut)
    expect_true(dir.exists(dOut), info="Test setup (output directory)")

    rSNVName <- tempfile()
    outFile <- .filter.utils.outputFilename(dOut, rSNVName)
    rSNVData <- c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t"),
        paste(c("TF2", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-02A-01W-0615-10",
            "8CD1CF38-DAD9-42EA-8B06-00955109F3D0"), collapse="\t"),
        paste(c("TF2", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0002-02A-01W-0615-10",
            "E6F5767F-E87B-4FC9-8880-CA70F434EEC9"), collapse="\t"),
        paste(c("HFE", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-02A-01W-0615-10",
            "F813F0B2-CE88-41FB-805E-40997E0E0309"), collapse="\t"),
        paste(c("HFE", "6", "127", "456", "+", "Missense_Mutation", "SNP",
            "A", "T", "G", "TCGA-A3-0002-02A-01W-0615-10", 
            "7719241D-B6C8-4B13-80F6-3047C8BBFE1F"), collapse="\t"))
    writeLines(rSNVData, con=rSNVName)
    expect_true(file.exists(rSNVName), info="Test setup (rSNV)")
    
    output <- filter.inconsistentCNA(rCNAName, c(rSNVName), dOut, filter)

    expect_equal(length(output), 2, info="2 genes filtered")
    expect_true('ZNF66' %in% output, info='ZNF66 has inconsistent sign')
    expect_true('TF2' %in% output, info='TF2 has inconsistent sign')

    fh <- file(outFile, open='r')
    filtered <- readLines(fh)
    close(fh)
    
    expect_equal(filtered, c(
        "#version 2.2",
        paste(c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", 
            "Strand", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode", "UUID"), collapse="\t"),
        paste(c("HFE", "6", "123", "456", "+", "Missense_Mutation", "SNP",
            "C", "T", "T",
            "TCGA-A3-0001-02A-01W-0615-10",
            "F813F0B2-CE88-41FB-805E-40997E0E0309"), collapse="\t"),
        paste(c("HFE", "6", "127", "456", "+", "Missense_Mutation", "SNP",
            "A", "T", "G", "TCGA-A3-0002-02A-01W-0615-10", 
            "7719241D-B6C8-4B13-80F6-3047C8BBFE1F"), collapse="\t")),
        info=paste(c("rSNV records with inconsistent sign removed", outFile)))

    genes <- readLines(con=filter)
    expect_equal(genes, output, info="output list same as filter file contents")
    
    file.remove(filter)
    file.remove(outFile)
    file.remove(rSNVName)
    file.remove(rCNAName)
    unlink(dOut)
})

