# testAnnotate.R
#
#
context("Gene graph annotation function")

# normal gH hashtable
# Unknown: "ENSP00000005257", "ENSP00000003084"
names = c("ENSP00000000233", "ENSP00000005340", "ENSP00000003100")
heats = c(0.3, 0.5, 0.2)
test_gH <- new.env(hash = TRUE)
for (i in 1:length(names)) {
    test_gH[[names[i]]] <- heats[i]
}

# build another hashtable with NULL heat value
normalNames = c("ENSP00000000233", "ENSP00000005340", "ENSP00000003100")
heatsWithNull = c(0.3, 0.5, NULL)
testgHNullHeat <- new.env(hash = TRUE)
for (i in 1:length(normalNames)) {
    testgHNullHeat[[normalNames[i]]] <- heatsWithNull[i]
}

test_gG <- importNet.STRING(fName = "dataString.txt",
                       net = "experimental",
                       dropUnmapped = FALSE,
                       silent = TRUE,
                       writeLog = FALSE)

test_that("annotateGraph rejects inappropriate objects", {
    # test null gH object
    null_gH <- NULL
    expect_error(annotateGraph(gH = null_gH, gG = test_gG, writeLog = FALSE),
                 "Object \"null_gH\" is NULL.")

    # test null gG object
    null_gG <- NULL
    expect_error(annotateGraph(gH = test_gH, gG = null_gG, writeLog = FALSE),
                 "Object \"null_gG\" is NULL.")

    # test gH object is missing
    expect_error(annotateGraph(gG = test_gG, writeLog = FALSE))

    # test gG object is missing
    expect_error(annotateGraph(gH = test_gH, writeLog = FALSE))

    # test writeLog is logical
    expect_error(annotateGraph(gH = test_gH, gG = test_gG, writeLog = "notlogical"),
                 "Error: writeLog must be of mode type and class logical")

    # test silent is logical
    expect_error(annotateGraph(gH = test_gH, gG = test_gG, silent = "notlogical"),
                 "Error: silent must be of mode type and class logical")
})

test_that("annotateGraph has added a new \"heat\" attribute", {
    test_AGG <- annotateGraph(gH = test_gH, gG = test_gG, writeLog = FALSE)
    # make sure it exists
    expect_true(!is.null(test_AGG))
    expect_true(!is.null(igraph::vertex_attr(test_AGG)$heat))
})

test_that("annotateGraph annotates graph vertices for which it has no score with 0", {
    test_AGG <- annotateGraph(gH = test_gH, gG = test_gG, writeLog = FALSE)
    expect_equal(igraph::vertex_attr(test_AGG)$heat[5], 0)
})

test_that("annotateGraph properly annotates the graph", {
    # check each heat attribute in the graph to make sure it's truthy
    test_AGG <- annotateGraph(gH = test_gH, gG = test_gG, writeLog = FALSE)
    expect_equal(igraph::vertex_attr(test_AGG)$heat, c(0.3, 0.5, 0.2, 0.0, 0.0))
})

test_that("annotateGraph returns a graph which is identical to the original gG,
          except for the heat attribute", {
    test_AGG <- annotateGraph(gH = test_gH, gG = test_gG, writeLog = FALSE)
    # test number of vertices
    expect_equal(igraph::vcount(test_gG), igraph::vcount(test_AGG))

    # test order of vertices
    expect_equal(igraph::vertex_attr(test_gG, "name"),
                 igraph::vertex_attr(test_AGG, "name"))

    # test number of edges
    expect_equal(igraph::ecount(test_gG), igraph::ecount(test_AGG))

    # test edge weights
    expect_equal(igraph::edge_attr(test_gG, "weight"),
                 igraph::edge_attr(test_AGG, "weight"))

    # test order of edges
    expect_equal(igraph::as_edgelist(test_gG), igraph::as_edgelist(test_AGG))
})

test_that("annotateGraph properly builds correct metadata for AGG", {
    test_AGG <- annotateGraph(gH = test_gH, gG = test_gG, writeLog = FALSE)
    expect_equal(attr(test_AGG, "type"), "AGG")
    expect_equal(attr(test_AGG, "version"), "1.0")
    expect_true(!is.null(attr(test_AGG, "UUID")))
})

# test_that("annotateGraph properly logs its event", {
#   # left out for now - needs a helper function since this will be
#   # tested often for other functions
# })

# think more about this one:
test_that("annotateGraph does not produce a broken AGG if it receives a broken input", {
    normal_AGG <- annotateGraph(gH = test_gH, gG = test_gG, writeLog = FALSE)

    # if there is a NULL value in heats
    test_null_heat_AGG <- annotateGraph(gH = testgHNullHeat, gG = test_gG,
                                        writeLog = FALSE)

    # that all heats are the same except for third heat value
    expect_equal(igraph::vertex_attr(test_null_heat_AGG)$heat[1:2],
                 igraph::vertex_attr(normal_AGG)$heat[1:2])
    expect_equal(igraph::vertex_attr(test_null_heat_AGG)$heat[4:5],
                 igraph::vertex_attr(normal_AGG)$heat[4:5])
})
