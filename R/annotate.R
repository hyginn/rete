# Annotate Graphs
#
#' Produce annotated gene graph where vertices have scores
#' attributed to genes.
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
#'   ## @examples ## \dontrun{
#'     annotateGraph(gH, gG)
#'   }
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

    # TODO: Need to check for environment variable, not really objects... how?
    # if (!exists(gH)) {
    #     stop(sprintf("Environment variable \"%s\" not yet set."), gHName)
    # }

    if (is.null(gG)) {
        stop(sprintf("Object \"%s\" is NULL.",
                     gGName))
    }

    # == EXTRACT ALL VERTICES =================================================
    allVertices <- igraph::vertex_attr(gG)$name

    # == CREATED ORDERED VECTOR OF HEAT ASSIGNMENTS ===========================
    numVerticesAnnotated <- 0
    heatAssignments <- c()
    for (i in 1:length(allVertices)) {
        if (!is.null(gH[[allVertices[i]]])) {
            heatAssignments[i] <- gH[[allVertices[i]]]
        } else {
            heatAssignments[i] <- 0
        }

        numVerticesAnnotated <- numVerticesAnnotated + 1

    }
    # == ATTACH CORRESPONDING VERTEX ANNOTATIONS TO AGG =========================
    AGG <- igraph::set_vertex_attr(gG, name = "heat", value = heatAssignments)

    # ==== SETUP METADATA ======================================================
    meta <- list(type = "AGG",
                 version = "1.0",
                 UUID = uuid::UUIDgenerate())

    if (!silent) {
        # cat()
    }

    if (writeLog) {
        logTitle <- "annotateGraph"

        # Compile function call record
        logCall <- character()
        logCall[1] <- "annotateGraph("
        logCall[2] <- sprintf("gH = \"%s\", ", as.character(gHName))
        logCall[3] <- sprintf("gG = \"%s\", ", as.character(gGName))
        logCall <- paste0(logCall, collapse = "")

        # Record progress information
        logNotes <- character()
        logNotes <- c(logNotes, sprintf("Annotated %s vertices", numVerticesAnnotated))

        # indicate output object name(s)
        logOutput = c("AGG")

        # send info to log file
        # logEvent(eventTitle = logTitle,
        #          eventCall = logCall,
        #          notes = logNotes,
        #          output = logOutput)
    }

    # ==== ATTACH METADATA =====================================================
    for (name in names(meta)) {
        attr(AGG, name) <- meta[[name]]
    }

    return(AGG)
}
