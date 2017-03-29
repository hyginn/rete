# importNet.STRING.R

#' Import network data from a STRING database file.
#'
#' \code{importNet.STRING} imports network edges from a STRING database file,
#' selects the highest confidence edges, maps ENSP IDs to HGNC gene symbols, and
#' returns a weighted, directed igraph graph object of a rete gG type.
#'
#' @section Selecting edges:
#'   STRING scores are p-values * 1000, rounded to
#'   integer. The function can retrieve  the highest scored edges according to
#'   three different cutoff type. Type "xN" (default: 10000) retrieves the xN
#'   highest scored edges. Type xQ (default 0.9) retrieves the edges with scores
#'   larger than the xQ quantile. Type "xS" (default 950) retrieves all edges
#'   with scores larger or equal to xS. If different values are requested, they
#'   are passed in the parameter val. To read all edges, cutoff type should be
#'   (the default) xN, val = Inf.
#' @section Networks:
#'   STRING "protein.links.detailed.v10.txt" files contain
#'   several protein networks: neighborhood, fusion, cooccurence, coexpression,
#'   experimental, database, textmining, and combined_score. However this
#'   function is not restricted to these types, but will read one network for
#'   which the column name is requested in the function's net parameter. This
#'   allows users to define their own column.
#' @section Tax ID:
#'   The taxID parameter is used as a sanity check on the file
#'   contents. Currently only the first data record protein IDs are checked for
#'   being prefixed with the tax ID. During processing all numbers and one
#'   period prefixed to the ENSP ID are removed.
#' @section Console and log:
#'   Progress is summarized to console and results are written to the log-file
#'   unless writeLog is FALSE.
#' @section The gG object:
#'   The function returns a rete gG object, a weighted,
#'   directed, simple igraph graph in which HGNC gene symbols are vertex names,
#'   and the edge attributes $weight hold the network scores. Metadata are
#'   stored as object attributs: $type: "gG"; $version: the gG
#'   object version; $UUID: the UUID which allows to retrieve detailed
#'   process information from the log file with \code{\link{findUUID}}.
#'
#' @param fName Filename of a STRING protein.links.detailed.v10.txt file.
#' @param net The requested network. This must be a string that exists in the
#'   header. The default is "combined_score".
#' @param cutoffType one of xN, xQ or xS (see Details).
#' @param val A number, quantile or score appropriate to the requested cutoff
#'   type.
#' @param taxID The NCBI tax ID prefix of the protein1 and protein2 IDs.
#'   Defaults to "9606" (homo sapiens)
#' @param dropUnmapped Controls whether to drop records in which at least one
#'  interactor could not be mapped to HGNC gene symbol. TRUE by default.
#'   by default.
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#' @return a weighted, directed, simple igraph graph which is a rete gG object.
#'
#' @family ImportNet.STRING, importNet.MultiNet, importNet.MITAB
#'
#'   ## @seealso \code{\link{fastMap}} fastMap() is used internally to map ENSP
#'   IDs to gene symbols.
#'
#'   ## @examples ## \dontrun{ ## importNet.STRING(IN, OUT) ## }
#' @export
importNet.STRING <- function(fName,
                             net = "combined_score",
                             cutoffType = "xN",
                             val,
                             taxID = "9606",
                             dropUnmapped = TRUE,
                             silent = FALSE,
                             writeLog = TRUE) {

    # ToDo: can we select the number of vertices ?
    # ToDo: complete writing progress to console
    # ToDo: Handle unmapped edges
    # ToDo: Remove duplicate edges and remove loops before subsetting,
    #         since this changes the number of edges above cutoff.
    # ToDo: write tests for that
    # ToDo: fix the fastMap() call
    # ToDo: stop on incomplete last line - suspect corrupt file

    #

    # ==== PARAMETERS ==========================================================

    cutoffTypes <- c("xS", "xQ", "xN")
    defaultValues <- c(xS = 950, xQ = 0.9, xN = 10000)

    if (missing(val)) {
        val <- defaultValues[cutoffType]
    }

    # ==== VALIDATIONS =========================================================

    # General parameter checks
    cR <- character()
    cR <- c(cR, .checkArgs(fName,        like = "FILE_E",   checkSize = TRUE))
    cR <- c(cR, .checkArgs(net,          like = "a",        checkSize = TRUE))
    cR <- c(cR, .checkArgs(cutoffType,   like = "a",        checkSize = TRUE))
    cR <- c(cR, .checkArgs(val,          like = numeric(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(taxID,        like = "a",        checkSize = TRUE))
    cR <- c(cR, .checkArgs(dropUnmapped, like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,     like = logical(1), checkSize = TRUE))

    if(length(cR) > 0) {
        stop(cR)
    }

    #  Validate the cutoffType parameter
    if (!cutoffType %in% cutoffTypes) {
        stop(paste0("Parameter error:", "\n",
                    "Valid cutoff types: ",
                    "\"", paste(cutoffTypes, collapse = "\" | \""), "\"", "\n",
                    "Found: \"", cutoffType, "\"", "\n"))
    }

    # Read header and one contents line
    tmp <- readLines(fName, n = 2)

    data   <- unlist(strsplit(tmp[2], " "))

    #  Validate theInteractor ID format
    patt <- paste(taxID, ".ENSP", sep = "")
    if (length(grep(patt, data[1:2])) != 2) {
        stop(paste0("ID error:", "\n",
                    "Expected format ", sprintf("<%s....>", patt), "\n",
                    "Found: <",paste(data[1:2], collapse = "> and <"), ">", "\n"))
    }

    # Parse for requested column. STRING data is " "-delimited.
    header <- unlist(strsplit(tmp[1], " "))
    iCol <- which(net == header)  # column that contains the requested net

    #  Validate that the requested network exists in data
    if (length(iCol) != 1) {
        stop(paste0("Request error:", "\n",
                    "Requested network type is ", sprintf("\"%s\"", net), "\n",
                    "Found network type(s): <",
                    paste(header[3:length(header)], collapse = ", "),
                    ">", "\n"))
    }

    # ==== READ DATA ===========================================================

    # Create column mask for readr::read_delim()
    readMask <- paste("cc",
                      paste(rep("_", length(header) - 2), collapse = ""),
                      sep = "")
    substr(readMask, iCol, iCol) <- "i"   # set target column to integer

    # start reading
    if (!silent) { cat("Reading", net, "from", fName, "...") }

    # Read file into data frame, use only selected column
    netDF <- readr::read_delim(file = fName,
                               delim = " ",
                               col_types = readMask,
                               n_max = Inf)

    colnames(netDF) <- c("a", "b", "weight")

    nIn <- nrow(netDF)
    if (!silent) { cat("done. (", nIn, " edges)", "\n") }


    # ==== SELECT EDGES ========================================================
    # Remove values below cutoff
    if (cutoffType == "xS") {
        # select all rows with weight >= val
        sel <- netDF$weight >= val
        nSel <- sum(sel)
    } else if (cutoffType == "xQ") {
        # select all rows with weight >= the val-quantile
        x <- stats::quantile(netDF$weight, probs = val)
        sel <- netDF$weight >= x
        nSel <- sum(sel)
    } else if (cutoffType == "xN") {
        # select the val highest scores, or all, whichever is fewer
        x <- order(netDF$weight, decreasing = TRUE)
        sel <- x[1:min(val, length(x))]
        nSel <- min(val, length(x))
    }

    netDF <- netDF[sel, ]  # Keep only selected edges

    # ==== MAP IDS TO GENE SYMBOLS =============================================

    # load fastMap hash table
    fastMapENSP <- readRDS(system.file("extdata",
                                       "fastMapENSP.rds",
                                       package = "rete"))

    # remove taxID prefix from STRING-format IDs
    netDF$a <- gsub("^[0-9]*\\.{0,1}", "", netDF$a)
    netDF$b <- gsub("^[0-9]*\\.{0,1}", "", netDF$b)

    # map ID to gene names
    mapA <- fastMap(netDF$a, fastMapENSP, type = "ENSP", dev = TRUE)
    mapB <- fastMap(netDF$b, fastMapENSP, type = "ENSP", dev = TRUE)

    if (dropUnmapped) {
        netDF$a <- mapA
        netDF$b <- mapB
        # select edges where one or the other gene symbol is NA
        sel <- is.na(netDF$a) | is.na(netDF$b)
        nDrop <- sum(sel)
        netDF <- netDF[-sel, ]  # Drop unmapped edges
    } else {
        # replace mapped IDs with their gene symbol
        netDF$a[! is.na(mapA)] <- mapA[! is.na(mapA)]
        netDF$b[! is.na(mapB)] <- mapB[! is.na(mapB)]
    }


    # ==== CREATE GRAPH OBJECT =================================================

    gG <- .df2gG(netDF)

    # ==== WRITE LOG ===========================================================

    if(writeLog) {

        myTitle <- "importNet.STRING"

        # Compile function call record
        myCall <- character()
        myCall[1] <- "importNet.STRING("
        myCall[2] <- sprintf("fname = \"%s\", ", fName)
        myCall[3] <- sprintf("net = \"%s\", ", net)
        myCall[4] <- sprintf("cutoffType = \"%s\", ", cutoffType)
        myCall[5] <- sprintf("val = %s, ", as.character(val))
        myCall[6] <- sprintf("taxID = \"%s\", ", taxID)
        myCall[7] <- sprintf("dropUnmapped = %s, ", as.character(dropUnmapped))
        myCall[8] <- sprintf("silent = %s, ", as.character(silent))
        myCall[9] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")

        # indicate input object name(s)
        #   (NA for this importNet.STRING())

        # Record progress information
        myNotes <- character()
        myNotes <- c(myNotes, sprintf("Read %s edges from file.", nIn))
        myNotes <- c(myNotes, sprintf("Selected %s edges via cutoff.", nSel))
        if (dropUnmapped) {
            myNotes <- c(myNotes, sprintf("Droppped %s unmapped edges.", nDrop))
        }
        myNotes <- c(myNotes, sprintf("gG object has %s vertices and %s edges.",
                                       igraph::gorder(gG), igraph::gsize(gG)))

        # indicate output object name(s)
        myOutput = c("gG")

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
#                input = character(),
                 notes = myNotes,
                 output = myOutput)
    }

    return(gG)
}


# [END]
