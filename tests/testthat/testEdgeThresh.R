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

# ((TODO: Is rete.EGGprototype only the barebones valid EGG, or does it contain testable
# data?))
# myEGG <- rete.EGGprototype
# Create test graph
myEGG <- igraph::random.graph.game(500, p.or.m = 1500, type = "gnm", directed = TRUE)
# Add edge weights
igraph::E(myEGG)$weight <- runif(n = igraph::ecount(myEGG), min = 0, max = 1)
# Size of the largest component in the test graph
nLargest <- max(igraph::components(myEGG)$csize)
# Largest component of myEGG
largestComponent <- igraph::decompose(myEGG, "weak")[[1]]
# Size, edge count, degrees and weights of the largest component
lcSize <- igraph::vcount(largestComponent)
lcEdgeCount <- igraph::ecount(largestComponent)
lcDegrees <- igraph::degree(largestComponent)
lcWeights <- igraph::E(largestComponent)$weight

# Permutated graph (|E| swaps)
permutedGraph <- .permuteGraph(myEGG, 1)
# Number of vertices in the test graph
nVertices <- igraph::vcount(myEGG)
# Number of edges in the test graph
nEdges <- igraph::ecount(myEGG)
# The degree of each vertex of the test graph
myDegrees <- igraph::degree(myEGG)
# ==== END SETUP AND PREPARE ===================================================


test_that("THRESH parameter errors are correctly handled", {
    # Try each parameter missing in turn (only EGG has no defaults)
    expect_error(THRESH(), info = "No Arguments Test")

    # Try parameters out of expected bounds
    expect_error(THRESH(myEGG, Lmax =  1), info = "Error: Lmax = 1")
    expect_error(THRESH(myEGG, Lmax =  0), info = "Error: Lmax = 0")
    expect_error(THRESH(myEGG, Lmax = -1), info = "Error: Lmax negative")

    expect_error(THRESH(myEGG, Lmin =  1), info = "Error: Lmin = 1 ")
    expect_error(THRESH(myEGG, Lmin =  0), info = "Error: Lmin = 0")
    expect_error(THRESH(myEGG, Lmin = -1), info = "Error: Lmin negative")

    expect_error(THRESH(myEGG, Lmax = 20, N =  0), info = "Error: N = 0")
    expect_error(THRESH(myEGG, Lmax = 20, N = -1), info = "Error: N negative")

    expect_error(THRESH(myEGG, Lmax = 20, Q =  0), info = "Error: Q = 0")
    expect_error(THRESH(myEGG, Lmax = 20, Q = -1), info = "Error: N negative")


    # Try NULL and length-zero parameters
    expect_error(THRESH(EGG = NULL,  Lmax = 20), info = "Error: EGG is NULL")
    expect_error(THRESH(EGG = myEGG, Lmax = 20, N = NULL), info = "Error: N is NULL")
    expect_error(THRESH(EGG = myEGG, Lmax = 20, Q = NULL), info = "Error: Q is NULL")
    expect_error(THRESH(EGG = myEGG, Lmax = NULL), info = "Error: Lmax is NULL")
    expect_error(THRESH(EGG = myEGG, Lmin = NULL), info = "Error: Lmin is NULL")

    # Commented out until .checkArgs is implemented
    # expect_error(THRESH(EGG = TRUE), info = "Error: EGG is boolean")
    # expect_error(THRESH(EGG = myEGG, Lmax = 20, N = TRUE), info = "Error: N is boolean")
    # expect_error(THRESH(EGG = myEGG, Lmax = 20, Q = TRUE), info = "Error: Q is boolean")
    # expect_error(THRESH(EGG = myEGG, Lmax = TRUE), info = "Error: Lmax is boolean")
    # expect_error(THRESH(EGG = myEGG, Lmin = TRUE), info = "Error: Lmin is boolean")
    #
    # expect_error(THRESH(EGG = "myEGG"), info = "Error: EGG is type character")
    # expect_error(THRESH(EGG =  myEGG , Lmax = "20"), info = "Error: Lmax is type character")
    # expect_error(THRESH(EGG =  myEGG , Lmin = "20"), info = "Error: Lmin is type character")
    # expect_error(THRESH(EGG =  myEGG , Q = "100"), info = "Error: Q is type character")
    # expect_error(THRESH(EGG =  myEGG , N = "100"), info = "Error: N is type character")

})

# ((Does .permuteGraph need parameter testing? or is it assumed that
# they are valid since the calling function has already checked them?))
test_that(".permuteGraph parameter errors are correctly handled", {

    expect_error(.permuteGraph(), info = "Error: No Arguments Given")
    expect_error(.permuteGraph(graph = myEGG), info = "Error: No Q Given")
    expect_error(.permuteGraph(Q = 100), info = "Error: No Graph Given")

    expect_error(.permuteGraph(graph = NULL,  Q = 100), info = "Error: Graph is NULL")
    expect_error(.permuteGraph(graph = myEGG, Q = NULL), info = "Error: Q is NULL")

    expect_error(.permuteGraph(graph = TRUE,  Q = 100), info = "Error: Graph is boolean")
    expect_error(.permuteGraph(graph = myEGG, Q = TRUE), info = "Error: Q is boolean")

    expect_error(.permuteGraph(graph ="myEGG",Q =  100 ), info = "Error: Graph is type character")
    expect_error(.permuteGraph(graph = myEGG, Q = "100"), info = "Error: Q is type character")

    expect_error(.permuteGraph(graph = myEGG, Q =  0), info = "Error: Q is 0")
    expect_error(.permuteGraph(graph = myEGG, Q = -1), info = "Error: Q is negative")

    expect_error(.permuteGraph(graph = 10,    Q = 100), info = "Error: Graph is type numeric")

})
# ((Should valid inputs different from rete.EGGprototype be used?))
# ((The testing below can get quite computationally expensive...))
test_that("a sane input to .permuteGraph gives an expected output", {

    # Check for self-loops and multiple edges
    expect_true(igraph::is.simple(permutedGraph), info = "Test: Permutation is simple")

    # Check that overall size, number of edges and node degrees of the
    # largest component has not changed
    expect_equal(igraph::vcount(permutedGraph), lcSize,
                 info = "Test: Permutation retains number of vertices")
    expect_equal(igraph::ecount(permutedGraph), lcEdgeCount,
                 info = "Test: Permutation retains number of edges in largest component")
    expect_equal(igraph::degree(permutedGraph), lcDegrees,
                 info = "Test: Permutation retains node degrees in largest component")

    # Check that weights are the same (Should they be checked edge by edge?)
    expect_true(all(igraph::E(myEGG)$weight %in% lcWeights),
                info = "Test: Permutation has the same edge weights")

})

test_that("a sane input to THRESH gives an expected output", {
    # Cover your code.

    ### Check THRESH with different valid input values

    gWithThresh <- THRESH(myEGG, Lmax = 20, N = 1, Q = 1)
    deltaVector <- igraph::graph.attributes(gWithThresh)$delta

    # Check that the input graph is not changed (by comparing all edges)
    # Get edge lists
    EL1 = igraph::get.edgelist(myEGG)
    EL2 = igraph::get.edgelist(gWithThresh)

    expect_true(all(duplicated(rbind(EL1, EL2))[(nrow(EL1) + 1):(2 * nrow(EL2))]),
                info = "Test: Input has the same edges as the output")

    # Check delta length, mode and range
    expect_equal(length(deltaVector), 19,
                 info = "Test: Delta attribute has expected length")
    expect_equal(mode(deltaVector),
                 "numeric", "Test: Delta attribute is numeric")
    expect_true(all(deltaVector[] >  0),
                info = "Test: Delta values are all positive")

    # (Are all delta values between 0 and 1?)
    # expect_true(all(deltaVector[] <= 1))

    # Check that weights are the same
    expect_true(all((igraph::E(myEGG)$weight %in% igraph::E(gWithThresh)$weight)))

    ###
    # Same as above, with Lmax of 10 (is this necessary?)
    gWithThresh <- THRESH(myEGG, Lmax = 10, N = 1, Q = 1)
    deltaVector <- igraph::graph.attributes(gWithThresh)$delta

    # Get edge lists
    EL1 = igraph::get.edgelist(myEGG)
    EL2 = igraph::get.edgelist(gWithThresh)

    expect_true(all(duplicated(rbind(EL1, EL2))[(nrow(EL1) + 1):(2 * nrow(EL2))]),
                info = "Test: Input has the same edges as the output")

    # Check delta length, mode and range
    expect_equal(length(deltaVector), 9)
    expect_equal(mode(deltaVector), "numeric")
    expect_true(all(deltaVector[] >  0))

    # (Are the delta values between 0 and 1?)
    # expect_true(all(deltaVector[] <= 1))

    # Check that weights are the same
    expect_true(all((igraph::E(myEGG)$weight %in% igraph::E(gWithThresh)$weight)))

})

test_that("a corrupt input to .permuteGraph does not lead to corrupted output", {
    # Expectations from .permuteGraph and THRESH functions are the same
    # in terms of handling corrupt inputs

    # ((TBD: Are non-simple graphs acceptable?))
    # invalidEGG1 <- igraph::rewire(myEGG, keeping_degseq(loops = TRUE))
    # invalidEGG2 <- igraph::rewire(myEGG, each_edge(
    #     loops = TRUE,
    #     multiple = TRUE,
    #     prob = 0.5
    # ))
    # expect_error(THRESH(EGG = invalidEGG1, Lmax = 20))
    # expect_error(THRESH(EGG = invalidEGG2, Lmax = 20))
})

test_that("a corrupt input to THRESH does not lead to corrupted output", {
    # Again: make absolutely sure be you never have an erroneous
    #        one-time execution of a loop

    # Test with invalid EGG graphs
    # ((Are non-simple graphs acceptable?))
    invalidEGG <- igraph::rewire(myEGG,
                                 igraph::each_edge(
                                     loops = TRUE,
                                     multiple = TRUE,
                                     prob = 0.99
                                 ))
    # Keep trying to create non-simple graphs
    while (!(igraph::any_multiple(invalidEGG)) ||
           !(any(igraph::which_loop(invalidEGG)))) {
        invalidEGG <- igraph::rewire(myEGG,
                                     igraph::each_edge(
                                         loops = TRUE,
                                         multiple = TRUE,
                                         prob = 0.99
                                     ))
    }
    expect_error(THRESH(
        EGG = invalidEGG,
        Lmax = 20,
        N = 1,
        Q = 1
    ),
    info = "Error: Non-Simple Input Graph")
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
rm(
    list = c(
        "largestComponent",
        "nLargest", "nEdges", "nVertices",
        "myEGG",
        "myDegrees",
        "permutedGraph",
        "lcDegrees",
        "lcSize",
        "lcEdgeCount",
        "lcWeights"))
# logName <- unlist(getOption("rete.logfile"))
# if (file.exists(logName)) { file.remove(logName)}
# options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
