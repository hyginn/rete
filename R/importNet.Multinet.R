# importNet.MULTINET.R

#' Import network data from BIOGRID, ENCODE data, KEGG, and SignaLink database 
#' (Multinet interactions).
#'
#' \code{importNet.MULTINET} imports network edges from a MULTINET database file,
#' 
#' @section Networks:
#'   MULTINET "protein.links.detailed.v10.txt" files contain
#'   several protein networks: neighborhood, fusion, cooccurence, coexpression,
#'   experimental, database, textmining, and combined_score. However this
#'   function is not restricted to these types, but will read one network for
#'   which the column name is requested in the function's net parameter. This
#'   allows users to define their own column.
#' @section The gG object:
#'   The function returns a rete gG object, a weighted,
#'   directed, simple igraph graph in which HGNC gene symbols are vertex names,
#'   and the edge attributes $weight hold the network scores. Graph attributes
#'   hold metadata; use igraph::graph_attr(gG) to return: $version: the gG
#'   object version; $logFile: the filename to which log information was
#'   written; $inFile: the input filename of MULTINET data; $call: the complete
#'   function call with expanded arguments; and $date: when the gG object was
#'   created.
#'
#' @param fName Filename of a MultiNet Multinet.interactions.network_presence.txt file.
#' @param net The requested network. This must be a MULTINET that exists in the
#'   header. The default is "PPI".
#' @param verbose Controls whether to print summary information after import. 
#'   Defaults to TRUE.
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
#'   ## @examples ## \dontrun{ ## importNet.MULTINET(IN, OUT) ## }
#' @export
importNet.MultiNet <- function(fName,
                             net = "PPI",
                             verbose = TRUE,
                             silent = FALSE,
                             writeLog = TRUE) {

    # ==== VALIDATIONS =========================================================

    # General parameter checks
    cR <- character()
    cR <- c(cR, .checkArgs(fName,        like = "FILE_E",   checkSize = TRUE))
    cR <- c(cR, .checkArgs(net,          like = "a",        checkSize = TRUE))
    cR <- c(cR, .checkArgs(verbose,      like = logical(1), checkSize = TRUE))
  
    if(length(cR) > 0) {
        stop(cR)
    }
    
    # Read header and one contents line
    tmp <- readLines(fName, n = 2)
    
    data   <- unlist(strsplit(tmp[2], " "))
    
    # Parse for requested column. STRING data is " "-delimited.
    header <- unlist(strsplit(tmp[1], " "))
    iCol <- which(net == header)  # column that contains the requested net
    
    #  Validate that requested network exists in data
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
    
    
    # ==== CREATE GRAPH OBJECT =================================================
    
    gG <- .df2gG(netDF)
    
    # ==== WRITE LOG ===========================================================
    
    if(writeLog) {
      
        myTitle <- "importNet.Multinet"
      
        # Compile function call record
        myCall <- character()
        myCall[1] <- "importNet.STRING("
        myCall[2] <- sprintf("fname = \"%s\", ", fName)
        myCall[3] <- sprintf("net = \"%s\", ", net)
        myCall[4] <- sprintf("verbose = %s, ", as.character(verbose))
        myCall[5] <- sprintf("silent = %s, ", as.character(silent))
        myCall[6] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")
      
        # indicate input object name(s)
        #   (NA for this importNet.Multinet())
      
        # Record progress information
        myNotes <- character()
        myNotes <- c(myNotes, sprintf("Read %s edges from file.", nIn))
        myNotes <- c(myNotes, sprintf("gG object has %s vertices and %s edges.",
                                       igraph::gorder(gG), igraph::gsize(gG)))
      
        # indicate output object name(s)
        myOutput = c("gG")
      
        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
                 notes = myNotes,
                 output = myOutput)
    }
    
    return(gG)
}


# [END]
