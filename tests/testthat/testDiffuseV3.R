# testDiffuseV3.R
#

context("Non-equilibrium heat diffusion")


# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()
# ==== END SETUP AND PREPARE ===================================================



fnAGG <- file.path(tempdir(), "AGG.rds")
writeLines("test", fnAGG)
fnEGG <- tempfile()

test_that("parameter errors are correctly handled", {
    # test fnAGG
    expect_error(
    diffuseV3(fnAGG = NULL,
              fnEGG = tempfile(),
              param = list(),
              silent = TRUE,
              writeLog = FALSE),
    '.checkArgs> "fnAGG" mode error')

    # test fnEGG
    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = NULL,
                  param = list(),
                  silent = TRUE,
                  writeLog = FALSE),
        '.checkArgs> "fnEGG" mode error')

    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = "no/such/dir/AGG.rds",
                  param = list(),
                  silent = TRUE,
                  writeLog = FALSE),
        '.checkArgs> "dirname\\(fnEGG\\)" error')

    # list of N, q, seed, useRCPP
    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = tempfile(),
                  param = list(N = numeric()),
                  silent = TRUE,
                  writeLog = FALSE),
        '.checkArgs> "N" length error')

    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = tempfile(),
                  param = list(q = numeric()),
                  silent = TRUE,
                  writeLog = FALSE),
        '.checkArgs> "q" length error')

    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = tempfile(),
                  param = list(seed = numeric()),
                  silent = TRUE,
                  writeLog = FALSE),
        '.checkArgs> "seed" length error')

    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = tempfile(),
                  param = list(useRcpp = numeric()),
                  silent = TRUE,
                  writeLog = FALSE),
        '.checkArgs> "useRcpp" length error')

    # Try parameters out of expected bounds
    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = tempfile(),
                  param = list(N = 0),
                  silent = TRUE,
                  writeLog = FALSE),
        'N must be greater then 0.')

    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = tempfile(),
                  param = list(q = -1),
                  silent = TRUE,
                  writeLog = FALSE),
        'q must be in the unit interval \\[0;1\\].')

    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = tempfile(),
                  param = list(q = 2),
                  silent = TRUE,
                  writeLog = FALSE),
        'q must be in the unit interval \\[0;1\\].')

    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = tempfile(),
                  param = list(seed = numeric()),
                  silent = TRUE,
                  writeLog = FALSE),
        '.checkArgs> "seed" length error')

    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = tempfile(),
                  param = list(useRcpp = numeric()),
                  silent = TRUE,
                  writeLog = FALSE),
        '.checkArgs> "useRcpp" mode error')

    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = tempfile(),
                  param = list(),
                  silent = logical(),
                  writeLog = FALSE),
        '.checkArgs> "silent" length error')

    expect_error(
        diffuseV3(fnAGG = fnAGG,
                  fnEGG = tempfile(),
                  param = list(),
                  silent = TRUE,
                  writeLog = logical()),
        '.checkArgs> "writeLog" length error')


    # Try NULL and length-zero parameters: be sure none lead to an erroneous
    #    one-time execution of any loop
    # Test NULL or zero length datastructures in .nonEqSim()
})


test_that("a sane input gives an expected output", {
    # Construct a small graph
    testDF <- data.frame(a = c("A", "A", "A", "B"),
                         b = c("B", "C", "D", "E"),
                         Weight = c(0.5, 1, 0, 1),
                         stringsAsFactors = FALSE)
    testAGG <- igraph::graph_from_data_frame(testDF, directed = TRUE)
    igraph::V(testAGG)$Score <- c(10, 0, 0, 0, 0)
    attr(testAGG, "type") <- "AGG"
    attr(testAGG, "version") <- "1.0"
    attr(testAGG, "UUID") <- uuid::UUIDgenerate()
    saveRDS(testAGG, fnAGG)

    diffuseV3(fnAGG,
              fnEGG,
              param = list(seed = 112358, useRcpp = FALSE),
              silent = TRUE,
              writeLog = FALSE)

    # Read result
    testEGG <- readRDS(fnEGG)
    expect_equal(names(igraph::V(testAGG)),            # Vertices equal
                 names(igraph::V(testEGG)))
    expect_equal(as.character(igraph::E(testAGG)),     # Edges equal
                 as.character(igraph::E(testEGG)))
    expect_equal(sum(igraph::V(testAGG)$Score),        # Total heat unchanged
                 sum(igraph::V(testEGG)$Score))
    expect_equal(sum(igraph::V(testEGG)$Score[c(2, 5)]), # Heat B + E
                 igraph::E(testEGG)$Influence[1])        # flux A->B
    expect_equal(igraph::V(testEGG)$Score[3],            # Heat C
                 igraph::E(testEGG)$Influence[2])        # flux A->C
    expect_equal(igraph::V(testEGG)$Score[4], 0)         # Heat D == 0
    # ToDo: Check code coverage
    # ToDo: Check log files

})

test_that("diffuseV3() doesn't leak output if silent = TRUE", {
    testF <- tempfile()
    capture.output(diffuseV3(fnAGG,
                             fnEGG,
                             param = list(seed = 112358, useRcpp = FALSE),
                             silent = TRUE,
                             writeLog = FALSE),
                   file = testF)
    expect_equal(length(readLines(testF)), 0)
    expect_true(file.remove(testF))
})


test_that("a corrupt input does not lead to corrupted output", {
    # ToDo
    # Try: - spurious characters in numeric columns ("N/A", ...).
    #      - extra tabs at line-end
    #      - comments before headers
    #      - truncated file (incomplete last line)
    # Again: make absolutely sure be you never have an erroneous
    #        one-time execution of a loop
})



# ==== BEGIN TEARDOWN AND RESTORE ==============================================

file.remove(fnAGG)
file.remove(fnEGG)

logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
