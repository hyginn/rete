# .importConservation.R

#' Utility functions from create and extracting UCSC Conservation Scores
.importConservation <- function(fName,
                                outName,
                                silent = FALSE) {

    # ==== PARAMETER TYPE VALIDATION ==========================================

    # General parameter checks
    cR <- character()
    cR <- c(cR, .checkArgs(fName,     like = "FILE_E",   checkSize = TRUE))
    cR <- c(cR, .checkArgs(outName,     like = "a",        checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,     like = logical(1),        checkSize = TRUE))

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

    currPos <- as.numeric(unlist(strsplit(header[startCol], "="))[2])

    # Determine the number of lines
    wcOut <- unlist(strsplit(system(paste("wc -l", fName), intern = TRUE), " +"))

    numLines <- as.numeric(wcOut[1])

    # Check if numLilnes is equal to
    if (numLines == 0) {
        errorMessage <- "Empty file."
        stop(errorMessage)
    }

    # Create a really big matrix
    chromosomeScores <- matrix(data = 0.0, ncol = 2, nrow = numLines)

    currIndex <- 1
    while (length(entry <- readLines(conservationFile, n = 1)) > 0) {
        if (!silent) {
            .pBar(currIndex, numLines)
        }
        # Read a header line
        # Scores are from 0.000 to 1.000
        # So if an entry was more than 5 characters in length, it's assumed to
        # be a header. Ten just ensures it is a header.
        if (nchar(entry) > 10) {
            header <- unlist(strsplit(entry, " "))
            currPos <- as.integer(unlist(strsplit(header[startCol], "="))[2])
        } else {
            num <- as.integer(as.numeric(entry) * 1000)
            chromosomeScores[currIndex, ] <- c(currPos, num)
            currIndex <- currIndex + 1
            currPos <- currPos + 1
        }
    }

    if (!silent) {
        .pBar(numLines, numLines)
    }

    close(conservationFile)

    # remove extra indices
    chromosomeScores <- chromosomeScores[1:(currIndex - 1), ]

    saveRDS(chromosomeScores, outName)
    rm(chromosomeScores)
}


# This is with the assumption that the geneStart, geneEnd are actually found in
# the chromosome data.
# This is also with an assumption that the there are no weird behaviour:
# start <- 1, end <- 100 but there is no data say from 20-30
# This can occur from multiple header lines in the UCSC chromosome score
# file.
.fetchConservationScores <- function(chromosomeScores,
                                   geneStart,
                                   geneEnd,
                                   cdsStarts,
                                   cdsEnds) {
    startIndex <- which(chromosomeScores[, 1] == geneStart)
    endIndex <- which(chromosomeScores[, 1] == geneEnd)
    extractedScores <- chromosomeScores[startIndex:endIndex, ]

    out <- matrix(nrow = 0, ncol = 2)
    # Assuming length(cdsStart) = length(cdsEnd)
    for (i in 1:length(cdsStarts)) {
        exonStart <- cdsStarts[i] - geneStart + 1
        exonEnd <- exonStart + (cdsEnds[i] - cdsStarts[i])
        exonScores <- extractedScores[exonStart:exonEnd, ]
        out <- rbind(out, exonScores)
    }
    return(out)
}
