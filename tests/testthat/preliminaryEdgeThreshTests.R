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
# ==== END SETUP AND PREPARE ===================================================


test_that("edgeThresh parameter errors are correctly handled", {
    # Try each parameter missing in turn (only EGG has no defaults)
    expect_error(edgeThresh())

    # Try parameters out of expected bounds
    # (Lmax, Q & N: 0 & negative, and Lmax = 1 & less than size of EGG's largest
    # component should result in error)


    # Try NULL and length-zero parameters: be sure none lead to an erroneous
    #   one-time execution of any loop: NULL, logical, character for: EGG,
    #   Lmax, N and Q

})

# ((Does .permuteGraph need parameter testing? or is it assumed that
# they are valid since the calling function has already checked them?))
# test_that(".permuteGraph parameter errors are correctly handled")
# Q = 0 and -1, NULL, logical or character should result in an error
# EGG = NULL, numeric, character or logical should result in an error
# Missing arguments should result in an error
test_that(".permuteGraph parameter errors are correctly handled") {



}
# ((Should valid inputs different from rete.EGGprototype be used?))
# ((The testing below can get quite computationally expensive...))
test_that("a sane input to .permuteGraph gives an expected output", {

    ### Test for Q = 1, 10, 20, 50, 100, 200
    # Set up will be as:
    myPermutation <- .permuteGraph(myEGG, Q = 1)

    # Check for self-loops and multiple edges
    expect_true(igraph::is.simple(myPermutation))

    # Check that the size of the largest connected component has not changed
    expect_equal(max(igraph::components(myPermutation)$csize), nLargest)

    # Check that overall size has not changed
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
    ### For each check, length deltaVector should be Lmax - 1, numeric, and
    ### between 0 and 1.
    ### Check that the input graph is not changed
    ### Test with N & Q = 1, 5, 10, Lmax 10 & 20
    ### ((Test with N = 10, Q = 1 -> is there a limit to how many permutations
    # there can be if only |E| edge swaps take place?))

    # First check will be like:

    gWithThresh <- edgeThresh(myEGG, Lmax = 20, N = 1, Q = 1)
    deltaVector <- igraph::graph.attributes(gWithThresh)$delta

    # Check that the input graph is not changed
    # ((what is the fastest way to check this?))
    expect_true(isomorphic(gWithThresh, myEGG))

    # Check delta length, mode and range
    expect_equal(length(deltaVector), 19)
    expect_equal(mode(deltaVector), "numeric")
    expect_true(all(deltaVector[] >  0))
    expect_true(all(deltaVector[] <= 1))



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
# Remove elements created during test setup

rm(
    list = c("myEGG",
             "nLargest",
             "permutedGraph",
             "nVertices",
             "nEdges"))

# logName <- unlist(getOption("rete.logfile"))
# if (file.exists(logName)) { file.remove(logName)}
# options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
