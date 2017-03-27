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
#'
#' @examples
#' \dontrun{combine(fNameSNV="SNV.rds", fNameCNA="CNA.rds", fgX="combined_data.rds", silent=TRUE, writeLog=FALSE)}
#'
#' @export
combineSNV_CNA <- function(fNameSNV, fNameCNA, fgX="gX.rds", silent=FALSE, writeLog=TRUE) {
    # parameter validation
    if (missing(fNameSNV) && missing(fNameCNA)) {
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
        stop("Output file has to be .rds extension!")
    }
    cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,     like = logical(1), checkSize = TRUE))
    if(length(cR) > 0) {
        stop(cR)
    }
    #
    # using hash package, e.g. extract keys later on
    hashTable <- hash()
    # read from input vector of SNV files
    for (i in 1:length(fNameSNV)) {
        SNVcurrent <- readRDS(fNameSNV[i])
        # get HGNC gene symbol row by row
        for (gene in rownames(SNVcurrent)) {
            position <- SNVcurrent[gene, "start"]
            vclass <- SNVcurrent[gene, "class"]
            # create key using <gene symbol>:SNV:<position>:<variant class>
            key <- sprintf("%s:SNV:%d:%s", gene, position, vclass)
            if (has.key(key, hashTable)) {
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
            # get CNA of one gene for each sample
            for (sample in colnames(CNAcurrent)) {
                tmpCNA <- CNAcurrent[gene, sample]
                # create key using <gene symbol:CNA:<as.character(round(<copy number>))>
                key <- sprintf("%s:CNA:%s", gene, as.character(round(tmpCNA)))
                if (has.key(key, hashTable)) {
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
    keys <- sort(keys(hashTable))
    # initialize the output dataframe
    out <- as.data.frame(matrix(ncol=4, nrow=length(keys)))
    colnames(out) <- c("sym", "type", "pos", "count")
    rownames(out) <- c(keys)

    for (k in keys) {
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
        logCall[2] <- sprintf("fNameSNV = \"%s\", ", fNameSNV)
        logCall[3] <- sprintf("fNameCNA = \"%s\", ", fNameCNA)
        logCall[4] <- sprintf("fgX = \"%s\", ", fgX)
        logCall[5] <- sprintf("silent = %s, ", as.character(silent))
        logCall[6] <- sprintf("writeLog = %s)", as.character(writeLog))
        logCall <- paste0(logCall, collapse = "")

        # Record progress information
        logNotes <- character()
        logNotes <- c(logNotes, sprintf("Combined %s samples in total.", length(hashTable)))

        # # send info to log file
        logEvent(eventTitle = logTitle,
                 eventCall = logCall,
                 notes = logNotes)
    }
}

# [END]
