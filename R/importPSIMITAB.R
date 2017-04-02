# importPSIMITAB.R

#' Import network data from a PSIMITAB database file.
#'
#' \code{importcore.psimitab} imports network edges from a PSIMITAB database file,
#' selects the highest confidence edges, maps UniProt IDs to HGNC gene symbols, and
#' returns a weighted, directed igraph graph object of a rete gG type.
#'
#' @section Selecting edges:
#'   PSIMITAB scores are p-values * 1000, rounded to
#'   integer. The function can retrieve  the highest scored edges according to
#'   three different cutoff type. Type "xN" (default: 10000) retrieves the xN
#'   highest scored edges. Type xQ (default 0.9) retrieves the edges with scores
#'   larger than the xQ quantile. Type "xS" (default 950) retrieves all edges
#'   with scores larger or equal to xS. If different values are requested, they
#'   are passed in the parameter val. To read all edges, cutoff type should be
#'   (the default) xN, val = Inf.
#' @section Networks:
#'   PSIMITAB "core.psimitab" files contain
#'   Unique gene identifiers for interactors. Alternative identifiers and aliases 
#'   for interactor genes. Interaction detection methods. Publication names. 
#'   Identifiers of publication in which the interaction has been shown. NCBI
#'   taxonomy IDs. Interaction types. Interaction identifiers. And finally confidence
#'   scores.
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
#' @param fName Filename of a PSIMITAB core.psimitab file.
#' @param cutoffType one of xN, xQ or xS (see Details).
#' @param val A number, quantile or score appropriate to the requested cutoff
#'   type.
#' @param taxID The NCBI tax ID prefix of the protein1 and protein2 IDs.
#'   Defaults to "9606" (homo sapiens)
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
importNet.MITAB <- function(fName,
                             cutoffType = "xN",
                             val,
                             experimentType = getOptions("rete.defaultPPI"),
                             taxID = "9606",
                             silent = FALSE,
                             writeLog = FALSE) {

  # ==== PARAMETERS =========================================================
  cutoffTytpes <- c("xS", "xQ", "xN")
  defaultValues <- c(xS = 950, xQ = 0.9, xN = 10000)

  if (missing(val)) {
    val <- defaultValues[cutoffType]
  }


  # ==== VALIDATIONS =======================================================

  cR <- character()
  cR <- c(cR, .checkArgs(fName,        like = "FILE_E",   checkSize = TRUE))
  cR <- c(cR, .checkArgs(cutoffType,   like = "a",        checkSize = TRUE))
  cR <- c(cR, .checkArgs(val,          like = numeric(1), checkSize = TRUE))
  cR <- c(cR, .checkArgs(taxID,        like = "a",        checkSize = TRUE))
  cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
  cR <- c(cR, .checkArgs(writeLog,        like = logical(1), checkSize = TRUE))

  if(length(cR) > 0) {
    stop(cR)
  }

  # Validate the cutoffType parameter
  if (!cutoffType %in% cutoffTypes) {
    stop(paste0("Parameter error:", "\n",
                "Valid cutoff types: ",
                "\"", paste(cutoffTypes, collapse = "\" | \""), "\"", "\n",
                "Found: \"", cutoffType, "\"", "\n"))
  }

  # Read one content line
  tmp <- readLines("core.psimitab", n = 1)
  data   <- unlist(strsplit(tmp, "\t"))

  if(sum(grepl("uniprotkb:", data[1:2])) != 2){
      stop(paste0("ID error:", "\n", "expected uniprot ID \n" paste("found" data[1:2])))
  }

  # ==== READ DATA ========================================================================

  readfile <- read.delim("core.psimitab", sep = "\t", header = FALSE)
  names(readfile)[1] <- "a"
  names(readfile)[2] <- "b"
  names(readfile)[15] <- "weight"
  scores <- readfile$weight

  # Create single score without | operator
  i = 1
  weightlist <- NULL
  while(i <= length(scores)){
    variable <- unlist(strsplit(as.character(scores[i]), "\\|"))
    if(variable[1] != "-" && variable[2] != "-"){
      value <- as.numeric(variable[1]) * as.numeric(variable[2])
      weightlist[i] <- value
    } else weightlist[i] <- 1
    i = i + 1
  }
  # Simplify data into only the needed coloumns
  readfile$weight <- weightlist
  a <- which(colnames(readfile) == "a")
  b <- which(colnames(readfile) == "b")
  c <- which(colnames(readfile) == "weight")
  netStr <- readfile[c(a,b,c)]

  if (cutoffType == "xS") {
    # select all rows with weight >= val
    sel <- netStr$weight >= val
    nSel <- sum(sel)
  } else if (cutoffType == "xQ"){
    x <- stats::quantile(netDF$weight, probs = val)
    sel <- netStr$weight <= val
    nSel <- sum(sel)
  } else if (cutoffType == "xN") {
    # select the val highest scores, or all, whichever is fewer
    x <- order(netStr$weight, decreasing = TRUE)
    sel <- x[1:min(val, length(x))]
    nSel <- min(val, length(x))
  }
  readfile <- readfile[sel, ] # Keep only selected genes

  fastMapUniProt <- readRDS("fastMapUniProt.rds")

  #Remove prefix from PSIMITAB format IDs
  netStr$a <- gsub(".*:", "", netStr$a)
  netStr$b <- gsub(".*:", "", netStr$b)

  # map ID to gene names
  mapA <- fastMap(netStr$a, fastMapUniProt, type = "UniProt", dev = TRUE)
  mapB <- fastMap(netStr$b, fastMapUniProt, type = "UniProt", dev = TRUE)

  # Drop unmapped genes
  netStr$a <- mapA
  netStr$b <- mapB
  # select edges where one or the other gene symbol is NA
  sel <- is.na(netStr$a) | is.na(netStr$b)
  nDrop <- sum(sel)
  netDF <- netStr[-sel, ]  # Drop unmapped edges


  gG <- .df2gG(netStr)

  # ==== WRITE LOG ===========================================================

  if(writeLog) {

    myTitle <- "importNet.STRING"

    # Compile function call record
    myCall <- character()
    myCall[1] <- "importNet.STRING("
    myCall[2] <- sprintf("fname = \"%s\", ", fName)
    myCall[3] <- sprintf("cutoffType = \"%s\", ", cutoffType)
    myCall[4] <- sprintf("val = %s, ", as.character(val))
    myCall[5] <- sprintf("taxID = \"%s\", ", taxID)
    myCall[6] <- sprintf("silent = %s, ", as.character(silent))
    myCall[7] <- sprintf("writeLog = %s)", as.character(writeLog))
    myCall <- paste0(myCall, collapse = "")


    # indicate input object name(s)
    #   (NA for this importNet.STRING())

    # Record progress information
    myNotes <- character()
    myNotes <- c(myNotes, sprintf("Read %s edges from file.", nIn))
    myNotes <- c(myNotes, sprintf("Selected %s edges via cutoff.", nSel))
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



