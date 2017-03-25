# importFilterHypermutators.R

#' Filters out hypermutators
#'
#' \code{importFilterHypermutators} identifies the number of SNV and CNA
#' mutations for each sample and removes the samples from rSNV and rCNA
#' datasets. Default threshold for removal is 400 mutations per sample.
#'
#' @param rSNVFileIn A vector of local file names of rSNV
#' @param rCNAFileIn A vector of local file names of rCNA
#' @param dOut Directory to store output, defaults to getwd().
#' @param xS Mutation threshold. 400 by default.
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#' @param writeDetailedLog Flag for extra details about log. TRUE by default.
#'
#' @examples
#' \dontrun{
#'     importFilterHypermutators(rSNVFileIn, rCNAFileIn, dOut, xS)
#' }
#'
#' @export
importFilterHypermutators <- function(rSNVFileIn = c(),
                                      rCNAFileIn = c(),
                                      dOut = getwd(),
                                      xS = 400,
                                      silent = FALSE,
                                      writeLog = TRUE,
                                      writeDetailedLog = TRUE) {

    ## VALIDATE PARAMS ########
    # make sure that files actually exist
    for (file in rSNVFileIn) {
        r <- .checkArgs(file, like = "FILE_E", checkSize = TRUE)
        if (length(r) > 0) {
            stop(r)
        }
    }

    for (file in rCNAFileIn) {
        r <- .checkArgs(file, like = "FILE_E", checkSize = TRUE)
        if (length(r) > 0) {
            stop(r)
        }
    }

    # make sure that output directory is valid
    r <- .checkArgs(dOut, like = "DIR", checkSize = TRUE)
    if(length(r) > 0) {
        stop(r)
    }

    r <- c(r, .checkArgs(xS, like = 1, checkSize = TRUE))
    if (length(r) > 0) {
        stop(r)
    }

    # set up hashes
    hashTable <- new.env(hash = TRUE)
    totalSamples <- 0

    ## COUNT MUTATIONS ########
    # for each file in rCNAFileIn
    for (file in rCNAFileIn) {
        # readRDS, count gene CNAs and add to hash for sample @ CNA and nMUT
        rCNA <- readRDS(file)
        for (sample in colnames(rCNA)) {
            totalSamples <- totalSamples + 1
            # if key is in hash, increment CNA and nMUT
            if (!is.null(hashTable$sample)) {
                prevCNACount <- hashTable$sample$CNA
                prevTotalCount <- hashTable$sample$total
                assign(sample, list(CNA = prevCNACount + 1, SNV = 0, total = prevTotalCount + 1), envir = hashTable)
            }
            # question - does each gene-CNA count as one mutation? I'm guessing yes?
            assign(sample, list(CNA = 1, SNV = 0, total = 1), envir = hashTable)
        }
    }

    # for each file in rSNVFileIn
    for (file in rSNVFileIn) {
        # readRDS, increment counter for SNV and nMUT
        rSNV <- readRDS(file)
        for (sample in rSNV$Tumor_Sample_Barcode) {
            # need to substitute dashes with period for consistency in files
            sample <-  gsub("-", ".", sample)
            totalSamples <- totalSamples + 1
            # if key is in hash, increment SNV and nMUT
            if (!is.null(hashTable$sample)) {
                prevCNACount <- hashTable$sample$CNA
                prevSNVCount <- hashTable$sample$SNV
                prevTotalCount <- hashTable$sample$total
                assign(sample, list(CNA = prevCNACount, SNV = prevSNVCount + 1, nMUT = prevTotalCount + 1), envir = hashTable)
            }
            assign(sample, list(CNA = 0, SNV = 1, total = 1), envir = hashTable)
        }
    }

    ## ASSESS MUTATIONS AND LOG STATISTICS #######

    # log global statistics:
    numSamplesBothSNVAndCNA <- 0
    numSamplesOnlyCNA <- 0
    numSamplesOnlySNV <- 0
    numSamplesNoChange <- 0
    numSamplesExceedThresh <- 0
    numRemovedSamples <- 0
    removedSamples <- c()

    for (sample in ls(hashTable)) {
        # num samples with only CNA
        if (hashTable[[sample]]$SNV == 0) {
            numSamplesOnlyCNA <- numSamplesOnlyCNA + 1
        }
        # num samples with only SNV
        else if (hashTable[[sample]]$CNA == 0) {
            numSamplesOnlySNV <- numSamplesOnlySNV + 1
        }
        # num samples with both SNV and CNA
        else if (hashTable[[sample]]$CNA > 0 && hashTable[[sample]]$SNV > 0) {
            numSamplesBothSNVAndCNA <- numSamplesBothSNVAndCNA + 1
        }
        # num samples that exceeded threshold and need removal
        else if (hashTable[[sample]]$total > xS) {
            numRemovedSamples <- numRemovedSamples + 1
            removedSamples <- c(removedSamples, sample)
        }
        # num samples with no change
        else {
            numSamplesNoChange <- numSamplesNoChange + 1
        }
    }

    ## PROCESS FILES AND UPDATE METADATA ##########

    # for each rCNAFileIn
    for (file in rCNAFileIn) {
        # for each sample, if hashTable$sample$total > xS, remove it
        # update metadata (getUUID(object))
        # save new CNAFile
    }


    # for each rSNVFileIn
    for (file in rSNVFileIn) {
        # for each sample (need to gsub("\.", "-", sample)), if hashTable$sample$total > xS, remove it
        # update metadata (getUUID(object))
        # save new SNVFile
    }

    ## LOGGING ###############

    if(writeLog) {

        logTitle <- "importFilterHypermutators"

        # Compile function call record
        logCall <- character()
        logCall[1] <- "importFilterHypermutators("
        logCall[2] <- sprintf("rSNVFileIn = \"%s\", ", rSNVFileIn)
        logCall[3] <- sprintf("rCNAFileIn = \"%s\", ", rCNAFileIn)
        logCall[4] <- sprintf("dOut = \"%s\", ", dOut)
        logCall[5] <- sprintf("xS = \"%s\", ", as.character(xS))
        logCall[6] <- sprintf("silent = %s, ", as.character(silent))
        logCall[7] <- sprintf("writeLog = %s)", as.character(writeLog))
        logCall <- paste0(logCall, collapse = "")

        # Record progress information
        logNotes <- character()
        logNotes <- c(logNotes, sprintf("Removed %s of %s samples", numRemovedSamples, totalSamples))

        if (writeDetailedLog) {
            logNotes <- c(logNotes, sprintf("%s of samples had both SNV and CNA", numSamplesBothSNVAndCNA/totalSamples))
            logNotes <- c(logNotes, sprintf("%s of samples had only CNA", numSamplesOnlyCNA/totalSamples))
            logNotes <- c(logNotes, sprintf("%s of samples had only SNV", numSamplesOnlySNV/totalSamples))
            logNotes <- c(logNotes, sprintf("%s of samples had no change", numSamplesNoChange/totalSamples))
            logNotes <- c(logNotes, sprintf("%s of samples exceeded threshold of %s", numSamplesExceedThresh/totalSamples, xS))

            # Accumulate all removed samples
            for (sample in removedSamples) {
                logNotes <- c(logNotes, sprintf("%s was removed", sample))
            }
        }

        # # send info to log file
        logEvent(eventTitle = logTitle,
                 eventCall = logCall,
                 notes = logNotes)
    }

}

# [END]
