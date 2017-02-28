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
