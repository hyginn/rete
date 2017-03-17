# testAnnotate.R
#
#
context("Gene graph annotation function")

# Unknown: "ENSP00000005257", "ENSP00000003084"
names = c("ENSP00000000233", "ENSP00000005340", "ENSP00000003100")
heats = c(0.3, 0.5, 0.2)
test_gH <- new.env(hash = TRUE)
for (i in 1:length(names)) {
    test_gH[[names[i]]] <- heats[i]
}

test_gG <- importNet.STRING(fName = "dataString.txt",
                       net = "experimental",
                       dropUnmapped = FALSE,
                       silent = TRUE,
                       writeLog = FALSE)

test_that("annotateGraph rejects inappropriate objects", {
    # test null gH object
    null_gH <- NULL
    expect_error(annotateGraph(gH = null_gH, gG = test_gG), "Object \"null_gH\" is NULL.")

    # test null gG object
    null_gG <- NULL
    expect_error(annotateGraph(gH = test_gH, gG = null_gG), "Object \"null_gG\" is NULL.")

    # test gH object is missing
    expect_error(annotateGraph(gG = test_gG))

    # test gG object is missing
    expect_error(annotateGraph(gH = test_gH))
})

test_that("annotateGraph has added a new \"heat\" attribute", {
    test_AGG <- annotateGraph(gH = test_gH, gG = test_gG)
    # make sure it exists
    expect_true(!is.null(test_AGG))
    expect_true(!is.null(igraph::vertex_attr(test_AGG)$heat))
})

test_that("annotateGraph annotates graph vertices for which it has no score with 0", {
    test_AGG <- annotateGraph(gH = test_gH, gG = test_gG)
    expect_equal(igraph::vertex_attr(test_AGG)$heat[5], 0)
})

test_that("annotateGraph properly annotates the graph", {
    # check each heat attribute in the graph to make sure it's truthy
    test_AGG <- annotateGraph(gH = test_gH, gG = test_gG)
    expect_equal(igraph::vertex_attr(test_AGG)$heat, c(0.3, 0.5, 0.2, 0.0, 0.0))
})

# think more about this one:
test_that("annotateGraph does not produce a broken AGG if it receives a broken input", {

})
