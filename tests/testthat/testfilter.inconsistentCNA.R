context("Filter CNA with inconsistent signs")

test_that("Filtering inconsistent CNA out of an rCNA RDS file", {
    outFile <- tempfile()
    rCNAName <- tempfile()
    filter <- tempfile()
    rCNAData <- data.frame(row.names = c("HFE", "FTO", "ZNF66", "TF2"),
        "TCGA-A3-0001-123" = c(3.0, -130.4, -1.5, 0.0),
        "TCGA-A3-0001-124" = c(2.9, -1.3, -0.5, 0.0),
        "TCGA-A3-0002-123" = c(1.0, -1.6, 3.0, 0.0),
        "TCGA-A3-0002-124" = c(1.0, -1.6, 3.0, 0.0))
    saveRDS(rCNAData, file=rCNAName)
    expect_true(file.exists(rCNAName), info="Test setup (rCNA)")
    
    output <- .filter.inconsistentCNA(rCNAName, outFile, filter)

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
})

