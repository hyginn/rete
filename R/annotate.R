# annotate.R

#' Produce annotated gene graph where vertices have scores
#' attributed to genes.
#'
#' \code{annotateGraph} takes as inputs a gene heat hashtable and a gene graph,
#' and annotates the graph with gene heats that have been calculated by the
#' scoring function. It then returns an AGG (annotated gene graph).
#'
#' @param gH List of genes and their associated score (mutational frequncy of
#'           mutSigCV value)
#' @param gG gene graph
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#' @return AGG object with annotated vertices
#'
#' @examples
#' \dontrun{
#'     annotateGraph(gH, gG)
#' }
#'
#' @export
annotateGraph <- function(gH, gG, silent = FALSE, writeLog = TRUE) {

    gHName <- deparse(substitute(gH))
    gGName <- deparse(substitute(gG))

    # check for params - that they are not missing
    if (is.null(gH)) {
        stop(sprintf("Object \"%s\" is NULL.",
                     gHName))
    }
    # Professor Steipe's .checkArgs()


    # TODO: Need to check for environment variable, not really objects... how?
    # if (!exists(gH)) {
    #     stop(sprintf("Environment variable \"%s\" not yet set."), gHName)
    # }

    if (is.null(gG)) {
        stop(sprintf("Object \"%s\" is NULL.",
                     gGName))
    }

    if (class(writeLog) != "logical") {
        stop("Error: writeLog must be of mode type and class logical")
    }

    if (class(silent) != "logical") {
        stop("Error: silent must be of mode type and class logical")
    }

    NL <- .PlatformLineBreak()

    # == EXTRACT ALL VERTICES =================================================
    if (!silent) { cat("Reading all vertices from ", gGName, NL) }

    allVertices <- igraph::vertex_attr(gG)$name

    if (!silent) { cat("Finished reading all vertices from ", gGName, NL) }

    # == CREATED ORDERED VECTOR OF HEAT ASSIGNMENTS ===========================
    if (!silent) { cat("Creating ordered vector of heat assignments", NL) }

    numVerticesAnnotated <- 0
    heatAssignments <- numeric(length(allVertices))

    for (i in 1:length(allVertices)) {
        if (!is.null(gH[[allVertices[i]]]) && !is.na(gH[[allVertices[i]]])) {
            heatAssignments[i] <- gH[[allVertices[i]]]
            numVerticesAnnotated <- numVerticesAnnotated + 1
        } else {
            heatAssignments[i] <- 0
        }
    }

    if (!silent) { cat("Finished creating ordered vector of heat assignments", NL) }

    # == ATTACH CORRESPONDING VERTEX ANNOTATIONS TO AGG =========================
    if (!silent) { cat("Attaching corresponding vertex annotations to AGG", NL) }

    AGG <- igraph::set_vertex_attr(gG, name = "heat", value = heatAssignments)

    if (!silent) { cat("Finished attaching corresponding vertex annotations to AGG", NL) }

    # ==== SETUP METADATA ======================================================
    meta <- list(type = "AGG",
                 version = "1.0",
                 UUID = uuid::UUIDgenerate())

    # ==== ATTACH METADATA =====================================================
    for (name in names(meta)) {
        attr(AGG, name) <- meta[[name]]
    }

    if (!silent) {
        cat()
    }

    if (writeLog) {
        logTitle <- "annotateGraph"

        # Compile function call record
        logCall <- character()
        logCall[1] <- "annotateGraph("
        logCall[2] <- sprintf("gH = \"%s\", ", gHName)
        logCall[3] <- sprintf("gG = \"%s\", ", gGName)
        logCall <- paste0(logCall, collapse = "")

        # Record progress information
        logNotes <- character()
        logNotes <- c(logNotes, sprintf("Annotated %s of %s vertices", numVerticesAnnotated, length(allVertices)))

        # indicate output object name(s)
        logOutput = c("AGG")

        # # send info to log file
        logEvent(eventTitle = logTitle,
                 eventCall = logCall,
                 notes = logNotes,
                 output = logOutput)
    }

    return(AGG)
}
