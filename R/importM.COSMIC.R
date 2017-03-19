# importM.COSMIC.R

#' Import mutation data from a COSMIC database.
#'
#' \code{importM.COSMIC}
#' Imports a tsv file containing CNV data from COSMIC. Imports a TSV file
#' mapping COSMIC gene IDs to HGNC gene symbols, required to match a CNV gene
#' symbol to an HGNC gene symbol (replaced by calls to fastMap()).
#' Returns an rCNA object in an RDS, where column names represents samples
#' (as TCGA-SampleID (note that the required information for a proper TCGA
#' mutation barcode is not available for COSMIC data, and that the sampleID
#' for this data is specific to COSMIC data)), row names are HGNC gene
#' symbols, and the values represent copy number.
#'
#' @section COSMIC CNV :
#' The format of the COSMIC CNV data is as follows.
#' (through sfpt, run get with this path:
#' /files/grch38/cosmic/v80/CosmicCompleteCNA.tsv.gz)
#'
#' [col number] Heading
#' [1] CNV_ID
#' [2] Gene ID
#' [4] Sample/Tumour ID
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
#' @section (Replaced with fastMap()) COSMIC HGNC:
#' The format of the COSMIC HGNC data is as follows.
#' (through sfpt, run get with this path:
#' /files/grch38/cosmic/v80/CosmicHGNC.tsv.gz)
#'
#' [col number] Heading
#' [1] COSMIC Gene ID
#' [2] COSMIC Gene Name
#' [3] Entrez ID
#' [4] HGNC ID
#' [5] Mutated? (does gene have coding mutations? (y/n))
#' [6] Cancer Cencus? (is gene in cancer cencus? (y/n))
#' [7] Expert Curated? (created by team of expert curators? (y/n))
#'
#' @param fNameCNV The path to a TSV file with CNV data from COSMIC.
#' #@param fNameHGNC The path to a TSV file with HGNC mapping data from COSMIC.
#' @param outFName The path to an output file for the rCNA RDS object. Only
#'   required if writeToFile == TRUE. Must be .rds file type.
#' @param writeToFile Controls whether to serialize rCNA RDS object.
#'   FALSE by default.
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#' @return Return rCNA object of COSMIC CNV data as RDS.
#'
#' @family importM.COSMIC, importM.TCGA, importM.GISTIC2, importM.ProjectGenie
#'
#'   ## @seealso \code{\link{fastMap}} fastMap() is used internally to
#'   map COSMIC IDs to gene symbols.
#'
#' @examples
#' \dontrun{
#' fNameCNV <- "./../test/COSMIC/testCosmicCNA.tsv"
#'
#' importM.COSMIC(fNameCNV,
#'                writeToFile = FALSE,
#'                silent = TRUE,
#'                writeLog = FALSE)
#' }
#' @export
importM.COSMIC <- function(fNameCNV,
                           #fNameHGNC,
                           outFName,
                           writeToFile = FALSE,
                           silent = FALSE,
                           writeLog = TRUE) {

    # fNameCNV <- "testCosmicCNA.tsv"
    # outFName <- "cosmicCNV.rds"
    # fNameHGNC <- "CosmicHGNC.tsv"

    # ==== PARAMETER TYPE VALIDATION ==========================================

    # Check that fNameCNV is a string.
    if(!is.character(fNameCNV)) {
        stop("fNameCNV is not a valid string.")
    }

    # # Check that fNameHGNC is a string.
    # if(!is.character(fNameHGNC)) {
    #     stop("fNameHGNC is not a valid string.")
    # }

    # Check that writeToFile is a logical.
    if (!is.logical(writeToFile)) {
        stop("writeToFile is not a valid logical.")
    }

    if (writeToFile) {
        # Make sure output file name is provided if writeToFile is TRUE.
        if (missing(outFName)) {
            stop("No output file specified in outFName.")
        }

        # Check that outFName is a string.
        if (!is.character(outFName)) {
            stop("outFName is not a valid string.")
        }

        # Check that outFName is of type .rds
        if (substr(outFName, nchar(outFName) - 3, nchar(outFName)) != ".rds") {
            stop("outFName is not of type .rds.")
        }
    } else {
        if (!missing(outFName)) {
            stop("outFName provided, but writeToFile is FALSE.")
        }
    }

    # Check that silent is a logical.
    if (!is.logical(silent)) {
        stop("silent is not a valid logicial.")
    }

    # Check that writeLog is a logical.
    if (!is.logical(writeLog)) {
        stop("writeLog is not a valid logicial.")
    }


    # ==== PARAMETER VALUE VALIDATION =========================================

    # Check that file provided at path fNameCNV exists
    if (!file.exists(fNameCNV)) {
        stop("Parameter error:\n   File specified at fNameCNV not found.")
    }

    # # Check that file provided at path fNameHGNC exists
    # if (!file.exists(fNameHGNC)) {
    #     stop("Parameter error:\n   File specified at fNameHGNC not found.")
    # }


    ## REPLACED BY CALLS TO fastMap() ##
    # # ==== READ COSMIC HGNC DATA AND CREATE COSMIC_ID->HGNC MAP================
    #
    # # Read header line from fNameHGNC file
    # tmpHeader <- readLines(fNameHGNC, n = 1)
    #
    # # Create character vector of col names from header string
    # header <- unlist(strsplit(tmpHeader[1], "\t"))
    #
    # # The required columns from Cosmic HGNC table
    # required <- c("COSMIC_ID", "HGNC_ID")
    #
    # # Find the column numbers in the header that match headings in required
    # iCol <- which(header %in% required)
    #
    # # Create a column mask for readr::read_delim()
    # readMask <- paste(paste(rep("_", length(header)),
    #                         collapse = ""),
    #                   sep = "")
    #
    # # Set target columns to characters
    # for (col in iCol) {
    #     substr(readMask, col, col) <- "c"
    # }
    #
    # # Read only the selected columns from Cosmic HGNC file
    # dataHGNC <- readr::read_delim(file = fNameHGNC,
    #                               delim = "\t",
    #                               col_types = readMask,
    #                               n_max = Inf)
    #
    # # Create a map of Cosmic IDs to HGNC IDs
    # HGNCmap <- dataHGNC[["HGNC_ID"]]
    # names(HGNCmap) <- dataHGNC[["COSMIC_ID"]]


    # ==== READ COSMIC CNV DATA ===============================================

    if (!silent) {
        cat("Reading CNV data from ", fNameCNV, " ...\n")
    }

    # Read header line from fNameCNV file
    tmpHeader <- readLines(fNameCNV, n = 1)

    # Create character vector of col names from header string
    header <- unlist(strsplit(tmpHeader[1], "\t"))

    # The required columns from Cosmic CNV table
    required <- c("ID_GENE", "ID_SAMPLE", "ID_TUMOUR", "TOTAL_CN")

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

    # Read only the selected columns from Cosmic CNV file
    dataCNV <- readr::read_delim(file = fNameCNV,
                                 delim = "\t",
                                 col_types = readMask,
                                 n_max = Inf)

    if (!silent) {
        cat("Done reading from ", fNameCNV, ".\n")
    }

    # ==== MAP COSMIC IDs TO HGNC GENE SYMBOLS ================================

    ## Used the following call to generate fastMapCosmic.rds:
    ## fastMapGenerate("hgnc_complete_set.txt", "symbol",
    ##                 "cosmic", type = "Cosmic",
    ##                 outputName = "../inst/extdata/fastMapCosmic.rds")

    # load fastMap hash table
    fastMapCosmic <- readRDS(system.file("extdata",
                                       "fastMapCosmic.rds",
                                       package = "rete"))

    # Validate that the type of the hash table is "Cosmic"
    if (!is.null(attributes(fastMapCosmic)$type)) {
        if ("Cosmic" != attributes(fastMapCosmic)$type) {
            stop(sprintf("Expected hash table type: Cosmic.\nSupplied hash table type: %s",
                         attributes(fastMapCosmic)$type))
        }
    } else {
        stop("Supplied hash table does not have a type attribute.")
    }


    # ==== CONSTRUCT rCNA OBJECT ==============================================

    if (!silent) {
        cat("Generating rCNA object using CNA data in ", fNameCNV, ".\n")
    }

    barcodes <- c()
    hgncSym <- c()

    for (i in 1:nrow(dataCNV)) {

        # Generate all barcodes for the mutations in dataCNV
        barcodes <- append(barcodes,
                           paste("TCGA-",
                                 dataCNV[i,]$ID_SAMPLE, "-",
                                 dataCNV[i,]$ID_TUMOUR,
                                 sep = ""))

        # Generate all HGNC symbols using fastMap()
        hgncSym <- append(hgncSym,
                          fastMap(dataCNV[i,]$ID_GENE,
                                  fastMapCosmic,
                                  type = "Cosmic",
                                  dev = TRUE))

    }

    # Create a matrix for storing CNV counts.
    #     nrow is equal to number of unique HGNC symbols
    #     ncol is equal to the number of unique barcodes.
    rCNA <- matrix(nrow = length(unique(hgncSym)),
                   ncol = length(unique(barcodes)))

    # Set the row and column names of rCNA to the
    # unique HGNC symbols and unique barcodes
    rownames(rCNA) <- unique(hgncSym)
    colnames(rCNA) <- unique(barcodes)

    # Populate rCNA matrix with CNA counts
    for (i in 1:nrow(dataCNV)) {
        rCNA[fastMap(dataCNV[i,]$ID_GENE,
                     fastMapCosmic,
                     type = "Cosmic",
                     dev = TRUE),
             paste("TCGA-",
                   dataCNV[i,]$ID_SAMPLE, "-",
                   dataCNV[i,]$ID_TUMOUR,
                   sep = "")] <- as.numeric(dataCNV[i,]$TOTAL_CN)
    }

    if (!silent) {
        cat("rCNA object successfully generated.\n")
    }


    # ==== WRITE TO LOG =======================================================

    if(writeLog) {

        myTitle <- "importM.COSMIC"

        # Compile function call record
        myCall <- character()
        myCall[1] <- "importM.COSMIC("
        myCall[2] <- sprintf("fnameCNV = \"%s\", ", fNameCNV)
        if (writeToFile) {
            myCall[3] <- sprintf("outFName = \"%s\", ", outFName)
        } else {
            myCall[3] <- ""
        }
        myCall[4] <- sprintf("writeToFile = %s, ", as.character(writeToFile))
        myCall[5] <- sprintf("silent = %s, ", as.character(silent))
        myCall[6] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")

        # Record progress information
        myNotes <- character()
        myNotes <- c(myNotes, sprintf("Read %s CNAs from file.", nrow(dataCNV)))
        if (writeToFile)
            myNotes <- c(myNotes, sprintf("Wrote rCNA object to %s.", outFName))

        # indicate output object name(s)
        # myOutput = c("rCNA")

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
                 notes = myNotes
                 # output = myOutput
        )
    }

    # Write to RDS file if specified
    if (writeToFile) {
        saveRDS(rCNA, outFName)
    }

    return(rCNA)

}

# [END]
