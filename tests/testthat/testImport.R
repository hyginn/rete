# testMultiply.R
#
#
context("import networks")

test_that("importNet.STRING reads STRING data", {
    gG <- importNet.STRING(fName = "dataString.txt",
                           net = "experimental",
                           silent = TRUE,
                           noLog = TRUE)
    expect_equal(igraph::vcount(gG), 5)
})

# [END]
