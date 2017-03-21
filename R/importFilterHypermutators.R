# importFilterHypermutators.R

#' Filters out hypermutators
#'
#' \code{importFilterHypermutators} identifies the number of SNV and CNA
#' mutations for each sample and removes the samples from rMUT and rCNA
#' datasets. Default threshold for removal is 400 mutations per sample.
#'
#' @param fNames A vector of local file names of rMUT and/or rCNA files
#' @param dOut Directory to store output, defaults to getwd().
#' @param xS Mutation threshold. 400 by default.
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#'
#' @examples
#' \dontrun{
#'     importFilterHypermutators(fNames, dOut, xS)
#' }
#'
#' @export
importFilterHypermutators <- function(fNames = c(),
                                      dOut = getwd(),
                                      xS = 400,
                                      silent = FALSE,
                                      writeLog = TRUE) {

    ## VALIDATE PARAMS ########
    # make sure that files actually exist
    for (file in fNames) {
        r <- .checkArgs(file, like = "FILE_E", checkSize = TRUE)
        if (length(r) > 0) {
            stop(sprintf("Local file name %s is invalid.", file))
        }
    }

    # make sure that output directory is valid
    r <- .checkArgs(fPath, like = "DIR", checkSize = TRUE)
    if(length(r) > 0) {
        stop(sprintf("Output directory %s is invalid.", dOut))
    }

    # make sure that xS is valid double
    r <- .checkArgs(xS, like = "double", checkSize = TRUE)
    if (length(r) > 0) {
        stop(xS)
    }


    # set up hash
    hashTable <- new.env(hash = TRUE)
    totalRemovedSamples <- 0
    totalSamples <- 0

    ## COUNT MUTATIONS ########
    # for each file in fNames,
        # determine if it is rMUT or rCNA
        # process rCNA files first
            # if rCNA: read file, count # of gene-CNAs for each sample, store
            # if rMUT: read file, if key not in hash then create
                    # increment counter of nMUT in sample's value

    ## ASSESS MUTATIONS #######
    # for each key in hash, set boolean thrsh <- TRUE if nMUT + nCNA > xS

    ## PROCESS FILES ##########

    # for each requested data file:
        # determine if rMUT or rCNA
            # if rCNA: read, remove columns where thrsh == TRUE, write file
            # if rMUT: read, do the same for new file

    ## LOGGING ###############

    if(writeLog) {

        logTitle <- "importFilterHypermutators"

        # Compile function call record
        logCall <- character()
        logCall[1] <- "importFilterHypermutators("
        logCall[2] <- sprintf("fNames = \"%s\", ", fNames)
        logCall[3] <- sprintf("dOut = \"%s\", ", dOut)
        logCall[4] <- sprintf("xS = \"%s\", ", as.character(xS))
        logCall[8] <- sprintf("silent = %s, ", as.character(silent))
        logCall[9] <- sprintf("writeLog = %s)", as.character(writeLog))
        logCall <- paste0(logCall, collapse = "")

        # Record progress information
        logNotes <- character()
        logNotes <- c(logNotes, sprintf("Removed %s of %s samples", totalRemovedSamples, totalSamples))

        # # send info to log file
        logEvent(eventTitle = logTitle,
                 eventCall = logCall,
                 notes = logNotes)
    }

    # returns nothing
}
