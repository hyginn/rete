# mapCoord.R

#' Populate a GeneData list with functions for mapping coordinates from chromosome to gene
#'
#' \code{mapCoord.R} populates map elements of a GeneData list with a
#' corresponding function closure that maps chromosomal coordinates to
#' gene coordinates.
#' The populated GeneData list is saved to an .rds file.
#'
#' @section Mapping function closure:
#' The mapping function that \code{mapCoord.R} creates and populates a GeneData
#' list which takes a chromosomal coordinate as input and returns the
#' corresponding gene coordinate if the chromosomal coordinate occurs in an
#' exon of the gene.
#' If the input is in an intron, the function returns 0.
#' If the input is out-of-bounds of the gene or invalid, the function
#' returns \code{NULL}.
#'
#' @param geneData the GeneData formatted resource list with at least strand,
#' start, end, exonStart, exonEnd elements populated
#' @param outFile File where the populated GeneData list is stored
#' (default: "geneData.rds")
#' @param dropExonBoundaries logical. Populated GeneData list has exonStart
#' and exonEnd elements dropped iff \code{TRUE}. (default: \code{TRUE})
#' @param silent logical. Console is suppressed iff silent is \code{TRUE}.
#' (default: \code{FALSE})
#' @param writeLog logical. Log file is written iff \code{TRUE}.
#' (default: \code{TRUE})
#'
#' @seealso \code{\link{importPDB}} Populates a GeneData list with 3D coordinates
#'
#' @export
mapCoord <- function(geneData,
                     outFile = "geneData.rds",
                     dropExonBoundaries = TRUE,
                     silent = FALSE,
                     writeLog = TRUE){

    #========================== Parameter validation ===========================
    cR <- character()
    cR <- c(cR, .checkArgs(geneData,           like = list()))
    cR <- c(cR, .checkArgs(outFile,            like = character()))
    cR <- c(cR, .checkArgs(dropExonBoundaries, like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,             like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,           like = logical(1), checkSize = TRUE))

    if(length(cR) > 0) {
        stop(cR)
    }

    # Parallelization for speed
    library(parallel)
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores, type = "PSOCK")
    on.exit(stopCluster(cl))

    # Parameter validation for geneData
    clusterExport(cl, c(".checkArgs", ".PlatformLineBreak"))

    parLapply(cl, geneData, function(x) {
        cRGeneData <- character()
        cRGeneData <- c(cRGeneData, .checkArgs(x$start,     like = integer(1),   checkSize = TRUE))
        cRGeneData <- c(cRGeneData, .checkArgs(x$end,       like = integer(1),   checkSize = TRUE))
        cRGeneData <- c(cRGeneData, .checkArgs(x$strand,    like = "+", checkSize = TRUE))
        cRGeneData <- c(cRGeneData, .checkArgs(x$exonStart, like = integer()))
        cRGeneData <- c(cRGeneData, .checkArgs(x$exonEnd,   like = integer()))
        if(length(cRGeneData) > 0) {
            stop(cRGeneData)
        }
    })

    #=============== Mapping function for a single gene model ==================
    makeMap <- function(ucscMod) {

        # Check if the number of exons listed in start and end match
        if(length(ucscMod$exonStart) != length(ucscMod$exonEnd)) {
            stop("Corrupt gene model")
        }

        # Prepare the environment (exon matrix, start, end, strand)
        exons <- cbind(ucscMod$exonStart, ucscMod$exonEnd)
        start <- ucscMod$start
        end <- ucscMod$end
        reverse <- FALSE
        if (ucscMod$strand == "-") {
            reverse <- TRUE
        }

        # Check if the exons and start stop points are validly ordered
        if(is.unsorted(c(start, as.vector(t(exons)), end))) {
            stop("Corrupt gene model")
        }

        # Generate cumulative sums
        exons <- cbind(exons, cumsum(exons[,2] - exons[,1] + 1))


        # Return the mapping function
        return(
            function(x) {
                # Check if input is numeric
                if (length(.checkArgs(x, like = integer(1), checkSize = TRUE)) > 0) {
                    stop("Non-numeric input")
                }

                # Check if out-of-range
                if (x < start | x > end) { return(NULL) }

                # Determine the exon our coordinate is located in
                iEx <- which(x >= exons[ , 1] & x <= exons[ , 2])

                # Check if intron
                if(length(iEx) != 1){ return(0) }

                # Exon bases to the left (toward 5', + strand gene start) of x
                total <- exons[iEx, 3] - (exons[iEx, 2] - x)

                if (reverse) {
                    #Exon bases to the right (toward 3', - strand gene start) of x
                    total <- exons[nrow(exons), 3] - total + 1
                }

                return(total)
            }
        )
    }

    clusterExport(cl, "makeMap", envir = environment())

    #====================== "Function manufacturing" ===========================
    geneData <- parLapply(cl, geneData, function(x) {
        x$map = makeMap(x)
        return(x)
    })

    if (dropExonBoundaries) {
        geneData <- parLapply(cl, geneData, function(x) {
            x$exonStart <- NULL
            x$exonEnd <- NULL
            return(x)
        })
    }

    saveRDS(geneData, file = outFile)

    #============================== Write log ==================================
    if(writeLog) {

        myTitle <- "mapCoord"

        # Compile function call record
        myCall <- character()
        myCall[1] <- "mapCoord("
        myCall[2] <- "geneData, "
        myCall[3] <- sprintf("outFile = \"%s\", ", outFile)
        myCall[7] <- sprintf("dropExonBoundaries = %s, ", as.character(dropExonBoundaries))
        myCall[8] <- sprintf("silent = %s, ", as.character(silent))
        myCall[9] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")

        # indicate input object name(s)
        myInput <- "geneData"

        # Record progress information
        myNotes <- sprintf("Generated %s mapping function closures", length(geneData))

        # indicate output object name(s)
        myOutput <- "geneData"

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
                 input = myInput,
                 notes = myNotes,
                 output = myOutput)
    }
}


# [END]
