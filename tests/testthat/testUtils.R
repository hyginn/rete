# testUtils.R
#
#
context("utility functions")
tmpFn <- tempfile()
writeLines("test", tmpFn)


test_that(".platformLineBreak() returns the correct linebreak", {
    if (.Platform$OS.type == "windows") {
        NL <- "\r\n"
    } else {
        NL <- "\n"
    }
    expect_equal(.PlatformLineBreak(), NL)
})


test_that(".checkArgs lets matching arguments pass", {
    # We are only counting the length of the error messsages, not contents
    expect_equal(length(.checkArgs(tempdir(), "DIR")), 0)
    expect_equal(length(.checkArgs(tmpFn, "FILE_E")), 0)
    expect_equal(length(.checkArgs(tmpFn, "FILE_W")), 0)

    uuid <- "e6c946b5-fd8c-4642-8c5c-a3736d5f5bab"
    expect_equal(length(.checkArgs(uuid, "UUID")), 0)
    uuid <- c(uuid, "{f6af5834-005d-40e0-afb7-dde9d438e489}")
    expect_equal(length(.checkArgs(uuid, "UUID")), 0)
    expect_equal(length(.checkArgs(toupper(uuid), "UUID")), 0)
    uuid <- "e6c946b5fd8c46428c5ca3736d5f5bab"
    expect_equal(length(.checkArgs(uuid, "UUID")), 0)

    expect_equal(length(.checkArgs(1, 1)), 0)
    expect_equal(length(.checkArgs(1L, 1:2)), 0)
    expect_equal(length(.checkArgs(as.matrix(1:2),
                                   as.matrix(3:4))), 0)
    expect_equal(length(.checkArgs(as.matrix(1:2),
                                   as.matrix(3:4), checkSize = TRUE)), 0)

    gG <- importNet.STRING(fName = "dataString.txt",
                           net = "experimental",
                           silent = TRUE,
                           writeLog = FALSE)
    gF <- gG
    expect_equal(length(.checkArgs(gF, gG)), 0)
})

test_that(".checkArgs finds argument errors", {
    # We are only counting the length of the error messsages, not contents
    expect_equal(length(.checkArgs("~/no/such/dir", "DIR")), 1)
    expect_equal(length(.checkArgs("noSuchFile", "FILE_E")), 1)
    expect_equal(length(.checkArgs("noSuchFile", "FILE_W")), 1)

    uuid <- "UUID"
    expect_equal(length(.checkArgs(uuid, "UUID")), 1)
    uuid <- c("e6c946b5-fd8c-4642-8c5c-a3736d5f5bab", "UUID")
    expect_equal(length(.checkArgs(uuid, "UUID")), 1)
    uuid <- c("UUID", "UUIE", "UUIF")
    expect_equal(length(.checkArgs(uuid, "UUID")), 3)
    uuid <- "00000000-0000-0000-0000-000000000000"    # reject NIL UUID
    expect_equal(length(.checkArgs(uuid, "UUID")), 1)

    expect_equal(length(.checkArgs(1, "1")), 3)  # mode, type and class error
    expect_equal(length(.checkArgs(1L, 1)), 2) # type and class error
    expect_equal(length(.checkArgs(1L, 1:2, checkSize = TRUE)), 1) # 1D
    expect_equal(length(.checkArgs(as.matrix(1:2),
                                   as.matrix(3:5), checkSize = TRUE)), 1)

    gG <- importNet.STRING(fName = "dataString.txt",
                           net = "experimental",
                           silent = TRUE,
                           writeLog = FALSE)
    gF <- gG
    igraph::graph_attr(gF)$newAttribute <- "test"  # add extra attribute
    expect_equal(length(.checkArgs(gF, gG)), 1) # number-of-attributes error

    gF <- gG  # revert
    igraph::graph_attr(gF)$version <- "noMatch"
    expect_equal(length(.checkArgs(gF, gG)), 1) # version error
})


test_that(".df2gG simplifies the graph", {
    netDF <- data.frame(a = c("crow", "duck", "crow", "crow"),
                        b = c("duck", "crow", "duck", "crow"),
                        weight = c(1, 2, 3, 4),
                        stringsAsFactors = FALSE)
    gG <- .df2gG(netDF)
    expect_equal(igraph::vcount(gG), 2)    # only duck and crow
    expect_equal(igraph::ecount(gG), 2)    # 1 duplicate, 1 self-edge removed
    expect_equal(igraph::edge_attr(gG)$weight, c(3, 2)) # correct weights
})


# [END]
