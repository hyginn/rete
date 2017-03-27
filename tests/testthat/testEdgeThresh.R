# edgeThresh.R
#
# Structure your tests according to the three principles expressed below:
#    test_that("parameter errors are correctly handled", { ... })
#    test_that("a sane input gives an expected output", { ... })
#    test_that("a corrupt input does not lead to corrupted output", { ...})


context("Graph permutation and delta calculation")


# ==== BEGIN SETUP AND PREPARE =================================================
# OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
# logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
# NL <- .PlatformLineBreak()
# if (file.exists(logName)) { file.remove(logName)}

# ((Is rete.EGGprototype only the barebones valid EGG, or does it contain testable
# data?))
myEGG <- rete.EGGprototype
nLargest <- max(igraph::components(EGG)$csize)
permutedGraph <- .permuteGraph(EGG)
nVertices <- igraph::vcount(EGG)
nEdges <- igraph::ecount(EGG)
myDegrees <- igraph::degree(EGG)
# ==== END SETUP AND PREPARE ===================================================


test_that("edgeThresh parameter errors are correctly handled", {
    # Try each parameter missing in turn (only EGG has no defaults)
    expect_error(edgeThresh())

    # Try parameters out of expected bounds
    expect_error(edgeThresh(myEGG, Lmax =  1))
    expect_error(edgeThresh(myEGG, Lmax =  0))
    expect_error(edgeThresh(myEGG, Lmax = -1))
    expect_error(edgeThresh(myEGG, Lmax = nLargst - 1))

    expect_error(edgeThresh(myEGG, Lmax = 20, N =  0))
    expect_error(edgeThresh(myEGG, Lmax = 20, N = -1))

    expect_error(edgeThresh(myEGG, Lmax = 20, Q =  0))
    expect_error(edgeThresh(myEGG, Lmax = 20, Q = -1))


    # Try NULL and length-zero parameters
    expect_error(edgeThresh(EGG = NULL,  Lmax = 20))
    expect_error(edgeThresh(EGG = myEGG, Lmax = 20, N = NULL))
    expect_error(edgeThresh(EGG = myEGG, Lmax = 20, Q = NULL))
    expect_error(edgeThresh(EGG = myEGG, Lmax = NULL))

    expect_error(edgeThresh(EGG = TRUE))
    expect_error(edgeThresh(EGG = myEGG, Lmax = 20, N = TRUE))
    expect_error(edgeThresh(EGG = myEGG, Lmax = 20, Q = TRUE))
    expect_error(edgeThresh(EGG = myEGG, Lmax = TRUE))

    expect_error(edgeThresh(EGG = "myEGG", Lmax =  20))
    expect_error(edgeThresh(EGG =  myEGG , Lmax = "20"))
    expect_error(edgeThresh(EGG =  myEGG , Lmax =  20, Q = "100"))
    expect_error(edgeThresh(EGG =  myEGG , Lmax =  20, N = "100"))

})

# ((Does .permuteGraph need parameter testing? or is it assumed that
# they are valid since the calling function has already checked them?))
test_that(".permuteGraph parameter errors are correctly handled", {

    expect_error(.permuteGraph())
    expect_error(.permuteGraph(graph = myEGG))
    expect_error(.permuteGraph(Q = 100))

    expect_error(.permuteGraph(graph = NULL,  Q = 100))
    expect_error(.permuteGraph(graph = myEGG, Q = NULL))

    expect_error(.permuteGraph(graph = TRUE,  Q = 100))
    expect_error(.permuteGraph(graph = myEGG, Q = TRUE))

    expect_error(.permuteGraph(graph ="myEGG",Q =  100 ))
    expect_error(.permuteGraph(graph = myEGG, Q = "100"))

    expect_error(.permuteGraph(graph = myEGG, Q =  0))
    expect_error(.permuteGraph(graph = myEGG, Q = -1))

    expect_error(.permuteGraph(graph = 10,    Q = 100))

})
# ((Should valid inputs different from rete.EGGprototype be used?))
# ((The testing below can get quite computationally expensive...))
test_that("a sane input to .permuteGraph gives an expected output", {

    ### Test for Q = 1, 10 as samples
    myPermutation <- .permuteGraph(myEGG, Q = 1)

    # Check for self-loops and multiple edges
    expect_true(igraph::is.simple(myPermutation))

    # Check that the size of the largest connected component has not changed
    expect_equal(max(igraph::components(myPermutation)$csize), nLargest)

    # Check that overall size, number of edges, node degrees has not changed
    expect_equal(igraph::vcount(myPermutation), nVertices)
    expect_equal(igraph::ecount(myPermutation), nEdges)
    expect_equal(igraph::degree(myPermutation), myDegrees)

    # Check that weights are the same
    expect_true(all((E(myEGG)$weight %in% E(myPermutation)$weight)))

    ###
    myPermutation <- .permuteGraph(myEGG, Q = 10)

    expect_true(igraph::is.simple(myPermutation))
    expect_equal(max(igraph::components(myPermutation)$csize), nLargest)
    expect_equal(igraph::vcount(myPermutation), nVertices)

})
# ((Is it possible to pass permuted graphs from the previous test section
# to the next one, if they pass? This would minimize a lot of computation))
test_that("a sane input to edgeThresh gives an expected output", {
    # Provide small input data, provide small output data.
    #   (Note that the output needs to be independently verifiably correct.
    #    Don't run an erroneous function, and then test against the erroneous
    #    output you received when you ran your function for the first time.)
    # Try this with important variations of parameters (reasonable, not
    #    combinatorially exhaustive).
    # Cover your code.

    # ((Is graph.attributes(gWithThresh)$delta the correct way of retrieving
    # delta values from graph attributes?)))
    ### Check edgeThresh with different valid input values

    gWithThresh <- edgeThresh(myEGG, Lmax = 20, N = 1, Q = 1)
    deltaVector <- igraph::graph.attributes(gWithThresh)$delta

    # Check that the input graph is not changed
    # ((what is the fastest way to check this?))
    expect_true(igraph::get.edgelist(myEGG) %in% get.edgelist(gWithThresh))

    # Check delta length, mode and range
    expect_equal(length(deltaVector), 19)
    expect_equal(mode(deltaVector), "numeric")
    expect_true(all(deltaVector[] >  0))
    expect_true(all(deltaVector[] <= 1))

    # Check that weights are the same
    expect_true(all((E(myEGG)$weight %in% E(gWithThresh)$weight)))

    ###
    gWithThresh <- edgeThresh(myEGG, Lmax = 10, N = 1, Q = 1)
    deltaVector <- igraph::graph.attributes(gWithThresh)$delta

    # Check that the input graph is not changed
    # ((what is the fastest way to check this?))
    expect_true(igraph::get.edgelist(myEGG) %in% get.edgelist(gWithThresh))

    # Check delta length, mode and range
    expect_equal(length(deltaVector), 9)
    expect_equal(mode(deltaVector), "numeric")
    expect_true(all(deltaVector[] >  0))
    expect_true(all(deltaVector[] <= 1))

    # Check that weights are the same
    expect_true(all((E(myEGG)$weight %in% E(gWithThresh)$weight)))

})

test_that("a corrupt input to .permuteGraph does not lead to corrupted output", {
    # Expectations from .permuteGraph and edgeThresh functions are the same
    # in terms of handling corrupt inputs
    invalidEGG1 <- igraph::rewire(myEGG, keeping_degseq(loops = TRUE))
    invalidEGG2 <- igraph::rewire(myEGG, each_edge(
        loops = TRUE,
        multiple = TRUE,
        prob = 0.5
    ))
    expect_error(edgeThresh(EGG = invalidEGG1, Lmax = 20))
    expect_error(edgeThresh(EGG = invalidEGG2, Lmax = 20))
})

test_that("a corrupt input to edgeThresh does not lead to corrupted output", {
    # Try: - spurious characters in numeric columns ("N/A", ...).
    #      - extra tabs at line-end
    #      - comments before headers
    #      - truncated file (incomplete last line)
    # Again: make absolutely sure be you never have an erroneous
    #        one-time execution of a loop

    # Test with invalid EGG graphs
    # ((How can one create an invalid EGG graph?))
    invalidEGG1 <- igraph::rewire(myEGG, keeping_degseq(loops = TRUE))
    invalidEGG2 <- igraph::rewire(myEGG, each_edge(
        loops = TRUE,
        multiple = TRUE,
        prob = 0.5
    ))
    expect_error(edgeThresh(EGG = invalidEGG1, Lmax = 20))
    expect_error(edgeThresh(EGG = invalidEGG2, Lmax = 20))
})



# ==== BEGIN TEARDOWN AND RESTORE ==============================================
rm(
    list = c(
        "gWithThresh",
        "deltaVector",
        "invalidEGG1",
        "invalidEGG2",
        "myPermutation",
        "myEGG",
        "permutedGraph"))
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
