# testMapCoord.R
#
#
context("GeneData mapping function generator")


# ==== BEGIN SETUP AND PREPARE =================================================

# Randomly generate a valid geneData structure of valid gene models
generateGeneData <- function() {
    geneData <- list()
    geneSymbols <- list("abcD", "efgH", "ijkL", "mnoP", "qrsT", "uvwX")

    # Convert a sorted even-sized list with a random number of random numbers
    # within an interval into a "gene model"
    for (i in 1:length(geneSymbols)) {
        lG <- 2 * sample(2:10, 1)
        preMod <- sample(1:10000, size=lG)
        preMod <- sort(preMod)
        plusmin <- c('+', '-')

        geneMod <- list(start=preMod[1], end=preMod[lG], strand=plusmin[sample(1:2, 1)])
        geneMod$exonStart <- c(preMod[2:(lG-1)])[c(TRUE, FALSE)]
        geneMod$exonEnd <- c(preMod[3:(lG)])[c(TRUE, FALSE)]

        geneData[[as.character(geneSymbols[i])]] <- geneMod
    }
    return(geneData)
}

# Build a new gene model based on running the mapFunctions on every coordinate
# If the mapping function was correct, the same gene model should be rebuilt
rebuildGeneModel <- function(mapFunction) {
    expGeneMod <- list()

    # Loop through all bases in the gene
    previous <- NULL
    counter <- 1
    for (i in 0:10001) {
        current <- mapFunction(i)
        if (is.null(previous) & !is.null(current)) {
            expGeneMod$start <- i
        } else if (!is.null(previous) & is.null(current)) {
            expGeneMod$end <- i - 1
        } else if (!is.null(previous) & !is.null(current)) {
            if (previous == 0 & current != 0) {
                expGeneMod$exonStart[counter] <- i
            } else if (previous != 0 & current == 0) {
                expGeneMod$exonEnd[counter] <- i - 1
                counter <- counter + 1
            }
        }
        previous <- current
    }

    return(expGeneMod)
}

# ==== END SETUP AND PREPARE ===================================================


test_that("parameter errors are correctly handled by mapCoord", {
    # Try setting each of the parameters to be NULL
    # Try setting "geneData" to a list with zero length lists

    expect_error(mapCoord())
    expect_error(mapCoord(NULL))
    expect_error(mapCoord(list(list(), list(), list())))

    gd <- generateGeneData()
    expect_error(mapCoord(gd, outFile = NULL))
    expect_error(mapCoord(gd, outFile = 42))
    expect_error(mapCoord(gd, dropExonBoundaries = NULL))
    expect_error(mapCoord(gd, dropExonBoundaries = "abcD"))
    expect_error(mapCoord(gd, silent = NULL))
    expect_error(mapCoord(gd, silent = list()))
    expect_error(mapCoord(gd, writeLog = NULL))
    expect_error(mapCoord(gd, writeLog = "42"))


})

test_that("a corrupt input does not lead to corrupted output for mapCoord", {
    # Generate geneData with generateGeneData and corrupt it.
    # Try corrupted "geneData":
    #   Exon starts and ends do not match
    #   Gene start and end do not match
    #   Types do not match

    wrongType <- generateGeneData()
    wrongType$abcD$exonEnd <- "Hello!"
    expect_error(mapCoord(wrongType))

    mismatchedExonsSingle <- generateGeneData()
    mismatchedExonsSingle$efgH$exonStart <- 1L
    expect_error(mapCoord(mismatchedExonsSingle))

    mismatchedExonsMulti <- generateGeneData()
    mismatchedExonsMulti$efgH$exonEnd <- 1:10
    expect_error(mapCoord(mismatchedExonsmismatchedExonsMulti))

    mismatchedStart <- generateGeneData()
    mismatchedStart$start <- 10001L
    expect_error(mapCoord(mismatchedStart))

    mismatchedEnd <- generateGeneData()
    mismatchedEnd$end <- 1L
    expect_error(mapCoord(mismatchedEnd))
})

test_that("parameter errors are correctly handled by the map functions produced by mapCoord", {
    # Generate geneData with generateGeneData and apply mapCoord to it
    # Try setting the map function parameters to NULL, non-ints

    gd <- generateGeneData()
    mapCoord(gd, writeLog = FALSE)
    gd <- readRDS("geneData.rds")
    for (name in names(gd)) {
        fx <- gd[[name]]$map
        expect_error(fx(NULL))
        expect_error(fx("42"))
        expect_error(fx(list()))
    }
})

test_that("the functions produced by mapCoord accurately map all coordinates", {
    # Generate geneData with generateGeneData
    # Check that mapCoord returns functions where the outputs of all the
    # coordinates can rebuild the correct gene model

    gd <- generateGeneData()
    mapCoord(gd, writeLog = FALSE)
    mappedGD <- readRDS("geneData.rds")
    for (name in names(gd)) {
        model <- gd[[name]]
        reModel <- rebuildGeneModel(mappedGD[[name]]$map)
        expect_equal(model$start, reModel$start)
        expect_equal(model$end, reModel$end)
        expect_equal(model$exonStart, as.vector(reModel$exonStart))
        expect_equal(model$exonEnd, as.vector(reModel$exonEnd))
    }
})

test_that("mapCoord drops exon boundaries iff dropExonBoundaries is TRUE", {
    # Generate geneData with generateGeneData
    # Test that dropExonBoundaries drops boundaries when TRUE, and default
    # Test that dropExonBoundaries does not drop boundaries when FALSE

    gd <- generateGeneData()

    mapCoord(gd, writeLog = FALSE)
    defaultDrop <- readRDS("geneData.rds")
    expect_equal(defaultDrop$abcD$exonStart, NULL)
    expect_equal(defaultDrop$abcD$exonEnd, NULL)

    mapCoord(gd, dropExonBoundaries = TRUE, writeLog = FALSE)
    trueDrop <- readRDS("geneData.rds")
    expect_equal(trueDrop$efgH$exonStart, NULL)
    expect_equal(trueDrop$efgH$exonEnd, NULL)

    mapCoord(gd, dropExonBoundaries = FALSE, writeLog = FALSE)
    falseDrop <- readRDS("geneData.rds")
    expect_false(is.null(falseDrop$uvwX$exonStart))
    expect_false(is.null(falseDrop$uvwX$exonEnd))
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
if (file.exists("geneData.rds")) { file.remove("geneData.rds") }
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
