# combineSNV_CNA.R

#' Combine mutation data from CNA and SNV
#'
#' \code{combineSNV_CNA} loads the vector of preprocessed dataset from filtered CNA or SNV files,
#'  creates hash keys corresponding to the input file type, combine the corresponding mutation count
#'  to each gene and save to a gX format file as output file.
#'
#' Details.
#' @section gX files:
#'  gX rds files contain data of SNV and CNA combined for each gene as one object, sorted per
#'  gene.
#'
#'
#' @param fNameCNA A vector of path of fully qualified file Names of SNV rds files.
#' @param fNameSNV A vector of path of fully qualifed filenames of CNA rds files.
#' @param fgX The fully qualified filename/path of a gX output file, by default is gX.rds.
#' @param silent Boolean option for writing combining process information to console, FALSE by default.
#' @param writeLog Boolean option for log results, TRUE by default.
#'
#' @return N/A. This function is normally invoked for its side effect of saving a gX file.
#'
#'
#' @examples
#' \dontrun{combineSNV_CNA(fNameSNV="SNV.rds", fNameCNA="CNA.rds", fgX="data.rds")}
#'
#' @export
combineSNV_CNA <- function(fNameSNV, fNameCNA, fgX="gX.rds", silent=FALSE, writeLog=TRUE) {
    # parameter validation
    if (missing(fNameSNV) && missing(fNameCNA)) {
        # ToDo: more descriptive error message
        stop("Input vectors are empty!")
    }
    cR <- character()
    for (fSNV in fNameSNV) {
        cR <- c(cR, .checkArgs(fSNV,        like = c("FILE_E")))
    }
    for (fCNA in fNameCNA) {
        cR <- c(cR, .checkArgs(fCNA,        like = c("FILE_E")))
    }
    if (grepl("\\.rds$", fgX) == FALSE) {
        # ToDo - need to discuss policy. I vote that this is not an error (bs)
        stop("Output file has to be .rds extension!")
    }
    cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,     like = logical(1), checkSize = TRUE))
    if(length(cR) > 0) {
        stop(cR)
    }
    #
    # using hash, e.g. extract keys later on
    hashTable <- hash::hash()
    # read from input vector of SNV files
    for (i in 1:length(fNameSNV)) {
        SNVcurrent <- readRDS(fNameSNV[i])
        # get HGNC gene symbol row by row
        for (gene in rownames(SNVcurrent)) {
            # ToDo - convert to i iterator and use progress bar.
            position <- SNVcurrent[gene, "start"]
            vclass <- SNVcurrent[gene, "class"]
            # create key using <gene symbol>:SNV:<position>:<variant class>
            key <- sprintf("%s:SNV:%d:%s", gene, position, vclass)
            if (hash::has.key(key, hashTable)) {
                hashTable[[key]] <- hashTable[[key]] + 1
            }
            else {
                hashTable[[key]] <- 1
            }
        }
    }
    # read from input vector of CNA files
    for (i in 1:length(fNameCNA)) {
        CNAcurrent <- readRDS(fNameCNA[i])
        # get HGNC gene symbols row by row
        for (gene in rownames(CNAcurrent)) {
            # ToDo - convert to i iterator and use progress bar.
            # get CNA of one gene for each sample
            for (sample in colnames(CNAcurrent)) {
                tmpCNA <- CNAcurrent[gene, sample]
                # create key using <gene symbol:CNA:<as.character(round(<copy number>))>
                # ToDo: this is a bit implicit. At this point we are changing the data
                #       format for CNAs from float to integer - and that's the value we
                #       actually record and use downstream. We need to make sure this is
                #       clearly expressed all over the docs.
                key <- sprintf("%s:CNA:%s", gene, as.character(round(tmpCNA)))
                if (hash::has.key(key, hashTable)) {
                    hashTable[[key]] <- hashTable[[key]] + 1
                }
                else {
                    hashTable[[key]] <- 1
                }
            }
        }
    }

    #
    # retrieve and sort all keys from hashTable
    keys <- sort(hash::keys(hashTable))
    # initialize the output dataframe
    out <- as.data.frame(matrix(ncol=4, nrow=length(keys)))
    colnames(out) <- c("sym", "type", "pos", "count")
    rownames(out) <- c(keys)

    for (k in keys) {
        # ToDo - convert to i iterator and use progress bar.
        # strplit each key to sym, type, pos and class and add to dataframe
        dataVector <- strsplit(k, ":")[[1]]

        sym <- dataVector[1]
        type <- dataVector[2]
        if (length(dataVector) == 3) {
            # create NA for $pos columns in CNAs
            pos <- NA
        }
        else {
            pos <- dataVector[3]
        }
        # retrieve value from hashTable using keys
        count <- hashTable[[k]]
        # add to dataframe
        out[k, "sym"] <- sym
        out[k, "type"] <- type
        out[k, "pos"] <- pos
        out[k, "count"] <- count
    }

    # update metadata
    # ToDo - add version.
    # ToDo - change type to "gX"
    # ToDo - CRITICAL add UUID. Saving the object with its old UUID
    #        is an intregration ERROR
    attr(out, "type") <- "combinedData"
    #attr(out, "UUID") <- getUUID(out)
    # save output file to local path

    saveRDS(out, file=fgX)

    #
    ## write log ###############

    if(writeLog) {

        logTitle <- "combineSNV_CNA"

        # Compile function call record
        logCall <- character()
        logCall[1] <- "combineSNV_CNA("
        # ToDo - this will not work as expected on multiple files
        logCall[2] <- sprintf("fNameSNV = \"%s\", ", fNameSNV)
        logCall[3] <- sprintf("fNameCNA = \"%s\", ", fNameCNA)
        logCall[4] <- sprintf("fgX = \"%s\", ", fgX)
        logCall[5] <- sprintf("silent = %s, ", as.character(silent))
        logCall[6] <- sprintf("writeLog = %s)", as.character(writeLog))
        logCall <- paste0(logCall, collapse = "")

        # Record progress information
        logNotes <- character()
        # ToDo: I would like to see separate summary values for CNAs and SNVs
        # ToDo: Review: if we save an object to file, and it has an UUID because
        #       it's an RDS, do we save the UUID to the logfile, even though this
        #       is not technically an "output" object and the UUID could change
        #       in the object.
        logNotes <- c(logNotes, sprintf("Combined %s samples in total.", length(hashTable)))

        # # send info to log file
        logEvent(eventTitle = logTitle,
                 eventCall = logCall,
                 notes = logNotes)
    }
}

# [END]
