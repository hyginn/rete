# importM.COSMIC.R

#' Import mutation data from a COSMIC database.
#'
#' \code{importM.COSMIC}
#' Returns a raw list of CNVs or NCVs from a COSMIC database.
#'
#' Details.
#' @section Requirements:
#' The file fname must be in TSV format,
#' and contain either CNV or NCV data from COSMIC.
#'
#' @param fName The path to a TSV file with CNV or NCV data from COSMIC.
#' @param type The type of data in fName - either "CNV" or "NCV".
#' @param range A vector of integers representing column ranges,
#'   specifying which columns from fName to return.
#' @param verbose Controls whether to print summary information after import.
#'   Defaults to TRUE.
#' @return Return selected columns from fName as a text file.
#'
#' @family importM.COSMIC, importM.TCGA, importM.GISTIC2, importM.ProjectGenie
#'
#' @examples
#' fName <- "./../test/COSMIC/CosmicNCV.tsv
#' type <- "NCV"
#' range <- c(1:2, 4)
#' verbose <- FALSE
#'
#' importM.COSMIC(fName, type, range, verbose)
#' @export
importM.COSMIC <- function(fName,
                           type,
                           range,
                           verbose = TRUE) {

    types <- c("CNV", "NCV") # fName must store either CNV or NCV data


    # ==== PARAMETER TYPE VALIDATION ==========================================

    # Check that fName is a string.
    if(!is.character(fName)) {
        stop("fName is not a valid string.")
    }

    # Check that type is a string.
    if(!is.character(type)) {
        stop("type is not a valid string.")
    }

    # If range is not provided, return all columns
    if (missing(range)) {
        noRange <- TRUE
    } else {
        noRange <- FALSE

        # Check that range is a vector.
        if(!is.vector(range)) {
            stop("range is not a valid vector.")
        }

        # Check that every index in range is a numeric.
        for (index in range) {
            if (!is.numeric(index)) {
                stop("range contains an invalid index.")
            }
        }

        # if range is empty, return all columns.
        if (length(range) == 0) {
            noRange <- TRUE
        }
    }

    if (!is.logical(verbose)) {
        stop("verbose is not a valid logicial.")
    }


    # ==== PARAMETER VALUE VALIDATION =========================================

    # Validate type parameter
    if (!type %in% types) {
        stop("Parameter error:\n   Valid file data types: \"",
             paste(types, collapse = "\" | \""),
             "\"\n   Found: \"",
             type,
             "\"\n")
    }

    # Check that file provided at path fName exists
    if (!file.exists(fName)) {
        stop("Parameter error:\n   File specified at fName not found.")
    }


    # ==== READ DATA ==========================================================

    # Read data from .tsv file into a data frame.
    cosmicData <- read.table(
        fName,
        sep="\t",
        header=TRUE,
        check.names = FALSE
    )


    # ==== FILTER COLUMNS =====================================================

    # Return all columns if range parameter was not provided.
    if (noRange) {
        cols <- c(1:ncol(cosmicData))
    # Otherwise, check that each column index in range is
    # greater than 0 and less than the number of columns in cosmicData.
    # If so, this is the range of columns needed from cosmicData
    } else {
        for (index in range) {
            if (index <= 0 || index > ncol(cosmicData)) {
                stop(paste("At least one index in range is out of bounds",
                           "(<= 0 or > number of columns in table",
                           "stored in fName file"))
            }
        }

        cols <- range
    }

    # Select columns from cosmicData based on indecies in cols
    selectedCols <- cosmicData[, cols]


    # ==== WRITE TO OUTPUT FILE ===============================================

    # Specify name of output text file based on type
    if (type == "CNV") {
        outFile <- "./CNV.txt"
    } else if (type == "NCV") {
        outFile <- "./NCV.txt"
    }

    # Write contents of selectedCols to output file.
    write.table(selectedCols,
                outFile,
                quote = FALSE,
                sep = "\t",
                row.names = FALSE)


    # ==== PRINT OUTCOME ======================================================
    if (verbose) {
        cat(paste("importM.COSMIC\n",
                    "Filtered columns in table from fName file.\n",
                    sprintf("Passed selected %d columns to %s",
                            length(cols),
                            outFile)))
    }
}
# [END]
