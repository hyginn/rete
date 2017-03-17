# testImport.R
#
#
context("import networks")

# ToDo test parameter handling
#      test correct logfile
#      test correct attributes

# pattern to grep for UUIDs
patt<-"[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}"

test_that("importNet.STRING produces gG from STRING data", {
    gG <- importNet.STRING(fName = "dataString.txt",
                           net = "experimental",
                           dropUnmapped = FALSE,
                           silent = TRUE,
                           writeLog = FALSE)
    expect_equal(igraph::vcount(gG), 5)         # correct number of vertices
    expect_equal(igraph::ecount(gG), 4)         # correct number of edges
    # metadata
    expect_equal(attr(gG, "type"), "gG")        # correct type
    expect_equal(attr(gG, "version"), "1.0")    # correct version
    expect_true(grepl(patt, attr(gG, "UUID")))  # valid UUID
})


# [END]
