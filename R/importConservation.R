# .importConservation.R

#' Generate a data file containing conservation score and chromosome position
#' pairs.
#'
#' \code{.importConservation} imports conservation scores from a UCSC chromosome
#' conservation score file. Generate a big matrix (col size 2, row size equal to
#' the number of rows of the input file) containing conservation score and
#' chromosome position pairs.
#'
#' @param fName The path to the UCSC chromosome conservation score file.
#' @param outName The path to an output file for the rCNA RDS object.
#' @param scoringType The type of the conservation score.
#' @param silent Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#' @param writeLog  Return rCNA object of COSMIC CNA data as RDS compressed data frame file.

.importConservation <- function(fName,
                                outName,
                                scoringType,
                                silent = FALSE,
                                writeLog = TRUE) {

    # ==== PARAMETER TYPE VALIDATION ==========================================

    # General parameter checks
    cR <- character()
    cR <- c(cR, .checkArgs(fName,     like = "FILE_E",   checkSize = TRUE))
    cR <- c(cR, .checkArgs(outName,     like = "a",        checkSize = TRUE))
    cR <- c(cR, .checkArgs(scoringType,     like = "a",        checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,     like = logical(1), checkSize = TRUE))

    if (length(cR) > 0) {
        stop(cR)
    }

    # ==== READ UCSC Conservation DATA ========================================

    if (!silent) {
        cat("Reading UCSC Conservation data from ", fName, " ...\n")
    }

    conservationFile <- file(fName, open = "r")

    # Read header file
    tmp <- readLines(conservationFile, n = 1)
    header <- unlist(strsplit(tmp, " "))

    chrCol <- grep("chr", header)
    startCol <- grep("start", header)

    if (length(chrCol) == 0 || length(startCol) == 0) {
        errorMessage <- "Invalid header."
        stop(errorMessage)
    }

    chromosomeSource <- unlist(strsplit(header[chrCol], "="))[2]
    currPos <- as.numeric(unlist(strsplit(header[startCol], "="))[2])
    currIndex <- 1


    # Determine the number of lines
    lineCountFile <- file(fName, open = "r")

    numLines <- 0

    while ((numLinesRead <- length(readLines(lineCountFile, n = 200000))) > 0) {
        numLines <- numLines + numLinesRead
    }

    # Check if numLilnes is equal to
    if (numLines == 0) {
        errorMessage <- "Empty file."
        stop(errorMessage)
    }

    close(lineCountFile)


    # Create really big matrix
    chromosomeScores <- matrix(data = 0.0, ncol = 2, nrow = numLines)

    while (length(entry <- readLines(conservationFile, n = 1)) > 0) {
        if (!silent) {
            .pBar(currIndex, numLines)
        }
        # Read a header line
        if (nchar(entry) > 10) {
            header <- unlist(strsplit(entry, " "))
            currPos <- as.integer(unlist(strsplit(header[startCol], "="))[2])
        } else {
            num <- as.integer(as.numeric(entry) * 1000)
            chromosomeScores[currIndex, 1] <- currPos
            chromosomeScores[currIndex, 2] <- num
            currIndex <- currIndex + 1
            currPos <- currPos + 1
        }
    }

    if (!silent) {
        .pBar(numLines - 1, numLines)
        .pBar(numLines, numLines)
    }

    close(conservationFile)

    # remove extra indices
    if (chromosomeScores[currIndex, 1] == 0) {
        currIndex <- currIndex - 1
    }
    chromosomeScores <- chromosomeScores[1:currIndex, ]

    # ==== SETUP METADATA ======================================================
    meta <- list(type = "ChromosomeConservationScore",
                 version = "1.0",
                 UUID = uuid::UUIDgenerate(),
                 chromosome = chromosomeSource,
                 scoreType = scoringType)


    # ==== ATTACH METADATA =====================================================

    for (name in names(meta)) {
        attr(chromosomeScores, name) <- meta[[name]]
    }

    # ==== WRITE TO LOG =======================================================

    if (writeLog) {

        myTitle <- ".importConservation"

        # Compile function call record
        myCall <- character()
        myCall[1] <- ".importConservaation("
        myCall[2] <- sprintf("fName = \"%s\", ", fName)
        myCall[3] <- sprintf("outName = \"%s\", ", outName)
        myCall[4] <- sprintf("scoringType = \"%s\", ", scoringType)
        myCall[5] <- sprintf("silent = %s, ", as.character(silent))
        myCall[6] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")

        # Record progress information
        myNotes <- character()
        myNotes <- c(myNotes, sprintf(
            "Read %s UCSC Conservation scores from file %s.",
            scoringType, fName))
        myNotes <- c(myNotes, sprintf(
            "Wrote Chromosome Conservation score to %s.", outName))

        # indicate output object name(s)
        myOutput = c("ChromosomeConservationScore")

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
                 notes = myNotes,
                 output = myOutput
        )
    }
    saveRDS(chromosomeScores, outName)
    rm(chromosomeScores)
}
