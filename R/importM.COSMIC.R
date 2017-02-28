# importM.COSMIC.R

#' Import mutation data from a COSMIC database.
#'
#' \code{importM.COSMIC} Returns a raw list of CNVs or NCVs from a COSMIC database.
#'
#' Details.
#' @section Requirements:
#' The file fname must be in TSV format, and contain either CNV or NCV data from COSMIC.
#'
#' @param fName The path to a TSV file with CNV or NCV data from COSMIC.
#' @param type The type of data in fName - either "CNV" or "NCV".
#' @param range A vector of string representing column ranges,
#'   specifying which columns from fName to return.
#' @param verbose Controls whether to print summary information after import. Defaults to TRUE.
#' @return Return selected columns from fName as a text file.
#'
#' @family importM.COSMIC, importM.TCGA, importM.GISTIC2, importM.ProjectGenie
#'
#' @examples
#' fName <- "./../test/COSMIC/CosmicNCV.tsv
#' type <- "NCV"
#' range <- c("1:2", "4")
#' verbose <- FALSE
#'
#' importM.COSMIC(fName, type, range, verbose)
#' @export
importM.COSMIC <- function(fName,
                           type,
                           range,
                           verbose = TRUE) {

    types <- c("CNV", "NCV") # fName must store either CNV or NCV data

    # ==== PARAMETER VALIDATION ===============================================

    # Validate type parameter
    if (!type %in% types) {
        stop("Parameter error:\n   Valid file data types: \"",
             paste(types, collapse = "\" | \""),
             "\"\n   Found: \"",
             type,
             "\"\n")
    }




}
