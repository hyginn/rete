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
#' @param outFName The path to an output file for the rCNA RDS object. Only
#'   required if writeToFile == TRUE. Must be .rds file type.
#' @param writeToFile Controls whether to serialize rCNA RDS object.
#'   FALSE by default.
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#' @return Return rCNA object of COSMIC CNA data as RDS.
#'
#' @family importM.COSMIC, importM.TCGA, importM.GISTIC2, importM.ProjectGenie
#'
#'   ## @seealso \code{\link{fastMap}} fastMap() is used internally to
#'   map COSMIC IDs to gene symbols.
#'
#' @examples
#' \dontrun{
#' fNameCNA <- "./../test/COSMIC/testCosmicCNA.tsv"
#'
#' importM.COSMIC(fNameCNA,
#'                writeToFile = FALSE,
#'                silent = TRUE,
#'                writeLog = FALSE)
#' }
#' @export
importM.COSMIC <- function(fNameCNA,
                           outFName,
                           writeToFile = FALSE,
                           silent = FALSE,
                           writeLog = TRUE) {


    # ==== PARAMETER TYPE VALIDATION ==========================================

    # General parameter checks
    cR <- character()
    cR <- c(cR, .checkArgs(fNameCNA,     like = "FILE_E",   checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeToFile,  like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,     like = logical(1), checkSize = TRUE))

    if(length(cR) > 0) {
        stop(cR)
    }

    if (writeToFile) {
        # Make sure output file name is provided if writeToFile is TRUE.
        if (missing(outFName)) {
            stop("No output file specified in outFName.")
        }

        # Check that outFName is a string.
        cR <- character()
        cR <- c(cR, .checkArgs(outFName,
                               like = "a",
                               checkSize = TRUE))

        if(length(cR) > 0) {
            stop(cR)
        }

    } else {
        if (!missing(outFName)) {
            stop("outFName provided, but writeToFile is FALSE.")
        }
    }

    # Check that file provided at path fNameCNA exists
    if (!file.exists(fNameCNA)) {
        stop("Parameter error:\n   File specified at fNameCNA not found.")
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

    # Set target columns to characters
    for (col in iCol) {
        substr(readMask, col, col) <- "c"
    }

    # Read only the selected columns from Cosmic CNA file
    dataCNA <- readr::read_delim(file = fNameCNA,
                                 delim = "\t",
                                 col_types = readMask,
                                 n_max = Inf)

    if (!silent) {
        cat("Done reading from ", fNameCNA, ".\n")
    }

    if (nrow(dataCNA) == 0) {
        stop("Input error: No CNAs found in file.")
    }

    # Convert all TOTAL_CN to numeric
    dataCNA$TOTAL_CN <- as.numeric(dataCNA$TOTAL_CN)

    # Subtract 2 from each CNA count
    if (!any(dataCNA$TOTAL_CN < 0)) {
        dataCNA$TOTAL_CN <-
            dataCNA$TOTAL_CN[dataCNA$TOTAL_CN >= 0] - 2
    }

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
        myCall[4] <- sprintf("writeToFile = %s, ", as.character(writeToFile))
        myCall[5] <- sprintf("silent = %s, ", as.character(silent))
        myCall[6] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")

        # Record progress information
        myNotes <- character()
        myNotes <- c(myNotes, sprintf("Read %s CNAs from file.", nrow(dataCNA)))
        if (writeToFile) {
            myNotes <- c(myNotes, sprintf("Wrote rCNA object to %s.", outFName))
        }

        # indicate output object name(s)
        myOutput = c("rCNA")

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
                 notes = myNotes,
                 output = myOutput
        )
    }

    # Write to RDS file if specified
    if (writeToFile) {
        saveRDS(rCNA, outFName)
    }

    return(rCNA)

}

# [END]
