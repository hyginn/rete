# testImport.R
#
#
context("import networks")

# ToDo test parameter handling
#      test correct logfile
#      test correct attributes

test_that("importNet.Multinet produces gG from Multinet data", {
    gG <- importNet.MULTINET(fName = "dataMultinet.txt",
                           net = "experimental",
                           silent = TRUE,
                           writeLog = FALSE)
    expect_equal(igraph::vcount(gG), 5)    # correct number of vertices
    expect_equal(igraph::ecount(gG), 4)    # correct number of edges
    expect_equal(attr(gG, "version"), "gG 1.0") # correct version
    expect_equal(attr(gG, "inFile"), "dataMultinet.txt") # metadata exists
})

test_that(".df2gG simplifies the graph", {
    netDF <- data.frame(a = c("crow", "duck", "crow", "crow"),
                        b = c("duck", "crow", "duck", "crow"),
                        weight = c(1, 2, 3, 4),
                        stringsAsFactors = FALSE)
    gG <- .df2gG(inFile = "dummy.txt", call = "dummy(arg = 1)")
    expect_equal(igraph::vcount(gG), 2)    # only duck and crow
    expect_equal(igraph::ecount(gG), 2)    # 1 duplicate, 1 self-edge removed
    expect_equal(igraph::edge_attr(gG)$weight, c(3, 2)) # correct weights
})



# [END]