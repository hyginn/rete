# importNet.STRING.R

#' Import network data from a STRING database file.
#'
#' \code{importNet.STRING} description.
#'
#' Details.
#' @section <title>: Additional explanation.
#'
#' @param <fName> <description>.
#' @param <net> <description>.
#' @param <cutoffType> <description>.
#' @param <val> <description>.
#' @param <taxID> <description>.
#' @param <silent> <description>.
#' @param <noLog> <description>.
#' @return <description>.
#'
#' @family ImportNet.STRING, importNet.MultiNet, importNet.MITAB
#'
#' ## @seealso \code{\link{<function>}} <describe related function>, ... .
#'
#' ## @examples
#' ## \dontrun{
#' ## importNet.STRING(IN, OUT)
#' ## }
#' @export
importNet.STRING <- function(fName,
                             net = "combined_score",
                             cutoffType = "xN",
                             val,
                             taxID = "9606",
                             silent = FALSE,
                             noLog = FALSE) {

    # ==== PARAMETERS ==========================================================

    cutoffTypes <- c("xS", "xQ", "xN")
    defaultValues <- c(xS = 950, xQ = 0.9, xN = 10000)

    if (missing(val)) {
        val <- defaultValues[cutoffType]
    }

    # Read header and 1 contents line
    tmp <- readLines(fName, n = 2)

    # Parse for requested column. STRING data is " "-delimited.
    header <- unlist(strsplit(tmp[1], " "))
    iCol <- which(net == header)  # column that contains the requested net

    data   <- unlist(strsplit(tmp[2], " "))


    # ==== VALIDATIONS =========================================================
    #    Validate cutoffType parameter
    if (!cutoffType %in% cutoffTypes) {
        stop("Parameter error:\n   Valid cutoff types: \"",
             paste(cutoffTypes, collapse = "\" | \""),
             "\"\n   Found: \"",
             cutoffType,
             "\"\n")
    }

    #    Validate Interactor ID format
    patt <- paste(taxID, ".ENSP", sep = "")
    if (length(grep(patt, data[1:2])) != 2) {
        stop("ID error:\n   Expected format ",
             sprintf("<%s....>", patt),
             "\n   Found: <",
             paste(data[1:2], collapse = "> and <"),
             ">\n")
    }

    #    Validate that requested network exists in data
    if (length(iCol) != 1) {
        stop("Request error:\n   Requested network type is ",
             sprintf("\"%s\"", net),
             "\n   Found network type(s): <",
             paste(header[3:length(header)], collapse = ", "),
             ">\n")
    }

    # ==== READ DATA ===========================================================

    # Create column mask for readr::read_delim()

    readMask <- paste("cc",
                      paste(rep("_", length(header) - 2), collapse = ""),
                      sep = "")
    substr(readMask, iCol, iCol) <- "i"   # set target column to integer

    # Read file into data frame, use only selected column
    netDF <- readr::read_delim(file = fName,
                               delim = " ",
                               col_types = readMask,
                               n_max = 10)

    colnames(netDF) <- c("a", "b", "weight")

    # Fix vertex names
    netDF$a <- gsub("^[0-9]*\\.{0,1}", "", netDF$a)
    netDF$b <- gsub("^[0-9]*\\.{0,1}", "", netDF$b)

    netDF$a <- fastMap(netDF$a, type = "ENSP")
    netDF$b <- fastMap(netDF$b, type = "ENSP")

    # Remove duplicates
    #   TODO

    # Remove values below cutoff
    #   TODO

    # Compile argument string
    fCall    <- "importNet.STRING("
    fCall[2] <- sprintf('fname = "%s", ', fName)
    fCall[2] <- sprintf("fname = \"%s\", ", fName)
    fCall[3] <- sprintf("net = \"%s\", ", net)
    fCall[4] <- sprintf("cutoffType = \"%s\", ", cutoffType)
    fCall[5] <- sprintf("val = %s, ", as.character(val))
    fCall[6] <- sprintf("taxID = \"%s\", ", taxID)
    fCall[7] <- sprintf("silent = %s, ", as.character(silent))
    fCall[8] <- sprintf("noLog = %s)", as.character(noLog))
    fCall <- paste(fCall, collapse = "")

    # ==== MAKE GRAPH ==========================================================
    gG <- .df2gG(fName, fCall, isDirected = TRUE)

    # ==== WRITE LOG ===========================================================
    if(! noLog) {
        logMessage    <- sprintf("%s: importNet.STRING()\n")
        logMessage[1] <- "    Returned gG object with"
        logMessage[2] <- sprintf("%d vertices and %d edges.\n",
                                 igraph::gorder(gG),
                                 igraph::gsize(gG))
        .appendToLog(paste(logMessage, collapse = ""))
    }

    return(gG)
}


# [END]
