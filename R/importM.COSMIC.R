# importM.COSMIC.R

#' Import mutation data from a COSMIC database.
#'
#' \code{importM.COSMIC}
#' Imports a tsv file containing CNA data from COSMIC. Imports a TSV file
#' mapping COSMIC gene IDs to HGNC gene symbols, required to match a CNA gene
#' symbol to an HGNC gene symbol (replaced by calls to fastMap()).
#' Returns an rCNA object in an RDS, where column names represents samples
#' (as TCGA-SampleID (note that the required information for a proper TCGA
#' mutation barcode is not available for COSMIC data, and that the sampleID
#' for this data is specific to COSMIC data)), row names are HGNC gene
#' symbols, and the values represent copy number.
#'
#' @section COSMIC CNA :
#' The format of the COSMIC CNA data is as follows.
#' (through sftp, run get with this path:
#' /files/grch38/cosmic/v80/CosmicCompleteCNA.tsv.gz)
#'
#' [col number] Heading
#' [1] CNA_ID
#' [2] Gene ID
#' [3] Gene Name
#' [4] Sample ID
#' [5] Tumour ID
#' [6] Primary Site
#' [7-9] Site Subtype 1-3
#' [10] Primary Histology
#' [11-13] Histology Subtype 1-3
#' [14] Sample Name
#' [15] Total Copy Number
#' [16] Minor Allele CN
#' [17] Mutation Type
#' [18] Study ID
#' [19] GRCh Coordinate System
#' [20] Genomic Coordinates of Variation
#'
#' @param fNameCNA The path to a TSV file with CNA data from COSMIC.
#' @param outFName The path to an output file for the rCNA RDS object.
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#' @return Return rCNA object of COSMIC CNA data as RDS compressed data frame file.
#'
#' @family importM.COSMIC, importM.TCGA, importM.GISTIC2, importM.ProjectGenie
#'
#'   ## @seealso \code{\link{fastMap}} fastMap() is used internally to
#'   map COSMIC IDs to gene symbols.
#'
#' @examples
#' \dontrun{
#' fNameCNA <- "./../test/COSMIC/testCosmicCNA.tsv"
#' outFName <- "cosmicCNA.rds"
#'
#' importM.COSMIC(fNameCNA,
#'                outFName,
#'                silent = TRUE,
#'                writeLog = FALSE)
#' }
#' @export
importM.COSMIC <- function(fNameCNA,
                           outFName,
                           silent = FALSE,
                           writeLog = TRUE) {


    # ==== PARAMETER TYPE VALIDATION ==========================================

    # General parameter checks
    cR <- character()
    cR <- c(cR, .checkArgs(fNameCNA,     like = "FILE_E",   checkSize = TRUE))
    cR <- c(cR, .checkArgs(outFName,     like = "a",        checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,     like = logical(1), checkSize = TRUE))

    if(length(cR) > 0) {
        stop(cR)
    }


    # ==== READ COSMIC CNA DATA ===============================================

    if (!silent) {
        cat("Reading CNA data from ", fNameCNA, " ...\n")
    }

    # Read header line from fNameCNA file
    tmpHeader <- readLines(fNameCNA, n = 1)

    # Create character vector of col names from header string
    header <- unlist(strsplit(tmpHeader[1], "\t"))

    # The required columns from Cosmic CNA table
    required <- c("gene_name", "ID_SAMPLE", "ID_TUMOUR", "TOTAL_CN")

    # Find the column numbers in the header that match headings in required
    iCol <- which(header %in% required)

    # Create a column mask for readr::read_delim()
    readMask <- paste(paste(rep("_", length(header)),
                            collapse = ""),
                      sep = "")

    # Set target columns to characters, except TOTAL_CN (read as integers)
    for (col in iCol) {
        if (header[col] == "gene_name") {
            substr(readMask, col, col) <- "c"
        } else {
            substr(readMask, col, col) <- "i"
        }
    }

    # Read only the selected columns from Cosmic CNA file
    dataCNA <- readr::read_delim(file = fNameCNA,
                                 delim = "\t",
                                 col_types = readMask,
                                 n_max = Inf)

    if (!silent) {
        cat("Done reading from ", fNameCNA, ".\n")
    }


    # ==== VALIDATE INPUT DATA ================================================

    # No gene_name should be any of these values
    illegalGeneNames <- c("NA", "N/A", "-")

    # Make sure the required columns are present in the input file.
    # Validated by checking the number of columns in dataCNA
    if (ncol(dataCNA) != 4) {
        stop("Input error: CNA file doesn't contain correct columns.")
    }

    # Make sure ID_SAMPLE, ID_TUMOUR, and TOTAL_CN are not "NA" after reading
    # as integer. This means there was a non-integer in these columns.
    # ToDo: Handle the case that there is a ligitimate missing value ...
    if (any(is.na(dataCNA$ID_SAMPLE)) ||
        any(is.na(dataCNA$ID_TUMOUR)) ||
        any(is.na(dataCNA$TOTAL_CN))) {
        stop(paste("Non-integer values found in ID_SAMPLE, " +
                   "ID_TUMOUR, and/or TOTAL_CN columns of input."))
    }

    # Make sure gene_name is not "N/A" or "NA" or "-" or NA.
    # ToDo: should probably use isHGNCsymbol() here...
    if (any(dataCNA$gene_name %in% illegalGeneNames) ||
        any(is.na(dataCNA$gene_name))) {
        stop("Illegal gene name found in input.")
    }

    # Make sure there is at least one CNA present in file
    if (nrow(dataCNA) == 0) {
        stop("Input error: No CNAs found in file.")
    }

    # Make sure all TOTAL_CN are non-negative. COSMIC CNA count data is not
    # expected to be negative. Must subtract 2 from each CNA count (done below)
    if (any(dataCNA$TOTAL_CN < 0)) {
        stop("Input error: CNA counts cannot be negative.")
    }


    # ==== FORMAT CNA DATA INPUT ==============================================

    # Subtract 2 from each CNA count. It is expected that TOTAL_CN counts
    # are not negative.
    # ToDo: Check back with author: conditional expression seems redundant with
    #       the check that was just made.
    dataCNA$TOTAL_CN <- dataCNA$TOTAL_CN[dataCNA$TOTAL_CN >= 0] - 2


    # ==== SETUP METADATA ======================================================

    meta <- list(type = "rCNA",
                 version = "1.0",
                 UUID = uuid::UUIDgenerate())


    # ==== CONSTRUCT rCNA OBJECT ==============================================

    if (!silent) {
        cat("Generating rCNA object using CNA data in ", fNameCNA, ".\n")
    }

    # Generate and get unique barcodes
    dataCNA$barcode <- paste("TCGA", dataCNA$ID_SAMPLE, dataCNA$ID_TUMOUR, sep = "-" )
    barcodes <- unique(dataCNA$barcode)

    # Get unique HGNC symbols
    hgncSym <- unique(dataCNA$gene_name)

    # Create a matrix for storing CNA counts.
    #     nrow is equal to number of unique HGNC symbols
    #     ncol is equal to the number of unique barcodes.
    rCNAMatrix <- (matrix(nrow = length(hgncSym),
                    ncol = length(barcodes)))

    rCNA <- data.frame(rCNAMatrix)

    # Set the row and column names of rCNA to the
    # unique HGNC symbols and unique barcodes
    rownames(rCNA) <- hgncSym
    colnames(rCNA) <- barcodes

    # Populate rCNA matrix with CNA counts
    for (i in 1:nrow(dataCNA)) {
        rCNA[dataCNA$gene_name[i], dataCNA$barcode[i]] <-
            dataCNA$TOTAL_CN[i]

        # Present a progress bar.
        if (!silent) {
            .pBar(i, nrow(dataCNA))
        }
    }


    # ==== ATTACH METADATA =====================================================

    for (name in names(meta)) {
        attr(rCNA, name) <- meta[[name]]
    }


    # ==== WRITE TO LOG =======================================================

    if(writeLog) {

        myTitle <- "importM.COSMIC"

        # Compile function call record
        myCall <- character()
        myCall[1] <- "importM.COSMIC("
        myCall[2] <- sprintf("fnameCNA = \"%s\", ", fNameCNA)
        myCall[3] <- sprintf("outFName = \"%s\", ", outFName)
        myCall[4] <- sprintf("silent = %s, ", as.character(silent))
        myCall[5] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")

        # Record progress information
        myNotes <- character()
        myNotes <- c(myNotes, sprintf("Read %s CNAs from file.", nrow(dataCNA)))
        myNotes <- c(myNotes, sprintf("Wrote rCNA object to %s.", outFName))

        # indicate output object name(s)
        myOutput = c("rCNA")

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
                 notes = myNotes,
                 output = myOutput
        )
    }

    # Write to RDS file
    saveRDS(rCNA, outFName)

}

# [END]
