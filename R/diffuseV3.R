# diffuseV3.R

#' Diffuse heat in an AGG object by numerical equilibration.
#'
#' \code{diffuseV3} Read an AGG object from an RDS compressed file. Diffuse
#'    a small portion q of excess heat across randomly chosen edges if the
#'    tail vertex of the directed edge has more "heat" than the head vertex.
#'    Repeat for N * |E| steps. Store the total amount of heat flowing across
#'    each edge as its "influence".
#'
#' @section Algorithm:
#'   Heat diffusion by numerical equilibration simulates the physical heat
#'   diffusion process by equilibrating portions of heat difference between
#'   adjacent vertices in directed, weighted graphs. The algorithm scales with
#'   O(|E|) and requires insignificant extra memory. It is thus well suited for
#'   graphs for which the calculation of heat diffusion from Markov chain
#'   stationary states via matrix inversion is too costly, as well as for
#'   weighted graphs. By default the seed for the random equilibration is taken
#'   from Sys.time(), but an explicit random seed can be passed in the parameter
#'   list. Heat diffusion is calculated in base R if useRcpp is FALSE, and uses
#'   Rcpp code otherwise (not yet implemented).
#'
#' @section Parameters:
#'   Typically the heat diffusion process is not run to equilibrium to emphasize
#'   the neighbourhood of the "hot" vertices in the graph. Large N and q values
#'   both speed up equilibration. Small q allow smooth equilibration without
#'   numerical artefacts. Thus q should be as small and N as large as
#'   computational resources allow. Setting a seed parameter allows to perform
#'   reproducible equilibrations for development purposes. If param is missing,
#'   the default is list(N = 100, q = 0.01, seed = as.numeric(Sys.time()),
#'   useRcpp = FALSE). If param is set, these default values will be used for
#'   any of the parameters that have not been specifed otherwise. ToDo: add
#'   benchmarks.
#'
#' @section Heat diffusion step:
#'   The amount of heat h that is transported
#'   between the two incident vertices of the directed edge e: (t->h), with h_t
#'   heat on the tail vertex and h_h heat on the head vertex, is: h <- q * W_e
#'   *(h_t - (h_t + h_h)/2) if h_t > h_h, else 0. q [0;1] is the global fraction
#'   of heat that is transported and governs the speed with which the graph will
#'   equilibrate to a state without heat differences; W_e [0;1] is the "thermal
#'   conductivity" of the specific edge. Appropriately scaled influence scores
#'   on edges enhance flow across these edges and throttle flow across others. h
#'   is subtracted from h_t and added to h_h thus the heat on both nodes
#'   asymptotically approaches the average from both sides. At every step the
#'   amount of heat that has flowed is added to the "Influence" attribute of the
#'   edge.
#'
#' @section Value:
#'   The resulting EGG graph is a igraph object, a directed,
#'   weighted graph with weights stored as a "Weight" edge attribute, the
#'   heat-flux as "Influence" edge-attributes, and the "heat" values as
#'   "Score" vertex attributes. Metadata is attached as attributes: "type":
#'   "EGG"; "version": "1.0"; and "UUID" a UUID.
#'
#' @section Log file:
#'   The filenames and UUID of the source AGG and the resulting EGG are written
#'   to the global log-file if writeLog is TRUE. As well, the function call
#'   and detailed parameters are written to the log-file. Thus the parameters
#'   with which an EGG was equilibrated can be found in the log-file by
#'   searching for its UUID.
#'
#' @param fnAGG Filename of an AGG object, saved as RDS.
#' @param fnEGG Filename of an EGG object, to be saved as RDS.
#' @param param list containing N, q, seed, and useRcpp parameters (see Details).
#' @param silent logical. Whether output will be to console will be suppressed.
#'                  Default FALSE.
#' @param writeLog logical. Whether an event entry for the log file is to be
#'                  written. Default TRUE.
#' @return N/A. This function is invoked for its side effect of writing a
#'   heat-equilibrated EGG graph to the ouput filename in RDS format.
#'
#'   ## @family ...
#'
#'   ## @seealso ...
#'
#' @examples
#' \dontrun{
#'    diffuseV3(IN, OUT)
#' }
#' @export
diffuseV3 <- function(fnAGG,
                        fnEGG,
                        param,
                        silent = FALSE,
                        writeLog = TRUE) {

    # ==== PARAMETERS ==========================================================

    if (missing(param)) {
        param <- list()
    }

    # use default values for any missing parameter
    if (is.null(param[["N"]])) {
        param[["N"]] <-100
    }
    if (is.null(param[["q"]])) {
        param[["q"]] <- 0.01
    }
    if (is.null(param[["seed"]])) {
        param[["seed"]] <- as.numeric(Sys.time())
    }
    if (is.null(param[["useRcpp"]])) {
        param[["useRcpp"]] <- TRUE
    }

    N       <- param[["N"]]
    q       <- param[["q"]]
    seed    <- param[["seed"]]
    useRcpp <- param[["useRcpp"]]

    # ==== VALIDATIONS =========================================================

    # General parameter checks
    cR <- character()
    cR <- c(cR, .checkArgs(fnAGG,          like = "FILE_E",  checkSize = TRUE))
    if (is.character(fnEGG)) {
        cR <- c(cR, .checkArgs(dirname(fnEGG), like = "DIR", checkSize = TRUE))
    } else {
        cR <- c(cR, .checkArgs(fnEGG,      like = "a",        checkSize = TRUE))
    }
    cR <- c(cR, .checkArgs(N,              like = 1,          checkSize = TRUE))
    cR <- c(cR, .checkArgs(q,              like = 1,          checkSize = TRUE))
    cR <- c(cR, .checkArgs(seed,           like = 1,          checkSize = TRUE))
    cR <- c(cR, .checkArgs(useRcpp,        like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,         like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,       like = logical(1), checkSize = TRUE))

    if(length(cR) > 0) {
        stop(cR)
    }

    # Range validations
    if (!N > 0) { stop("N must be greater then 0.")}
    if (!q >= 0 || !q <= 1 ) {
        stop("q must be in the unit interval [0;1].")
    }

    # ==== READ AGG ===========================================================

    # Read AGG (a weighted directed graph) as an RDS object.
    AGG <- readRDS(fnAGG)

    # Validate the graph
    if (!igraph::is_directed(AGG)) {
        stop(sprintf("The input graph read from % is not a directed graph.",
                     fnAGG))
    }

    # ==== SETUP DATSTRUCTURES =================================================
    # Setup a suitable data structure for edges and vertices.

    tmp <- igraph::ends(AGG, igraph::E(AGG), names = FALSE)
    eTable <- data.frame(a = tmp[ , 1],
                         b = tmp[ , 2],
                         W = igraph::E(AGG)$Weight,
                         flux = numeric(nrow(tmp)),
                         stringsAsFactors = FALSE)
    rm(tmp)
    vTable <- data.frame(heat = igraph::V(AGG)$Score)

    # Validate the weights
    if (!all(eTable$W >= 0 & eTable$W <= 1)) {
        stop(sprintf("One or more E(AGG)$Weight values %s",
                     "are not in the unit interval."))
    }

    # ==== RUN SIMULATION =================================================

    if (!silent) {
        cat(sprintf("Beginning equilibration for N * |E| = %d * %d %s.\n",
                    N,
                    nrow(eTable),
                    "steps"))
    }

    # Call appropriate heat equilibration function depending on "useRcpp".
    if (useRcpp) {
        stop("Panic: Rcpp code not yet implemented.")
    } else {

        EGGdata <- .nonEqSim(vTable,
                             eTable,
                             N,
                             q,
                             seed,
                             silent)

    }

    rm(eTable)
    rm(vTable)

    # reassemble function output to EGG (igraph object)

    EGG <- igraph::set_edge_attr(AGG, "Influence", value=EGGdata$eTable$flux)
    EGG <- igraph::set_vertex_attr(EGG, "Score",   value=EGGdata$vTable$heat)

    # attach metadata to EGG object
    attr(EGG, "type")    <- "EGG"
    attr(EGG, "version") <- "1.0"
    attr(EGG, "UUID")    <- uuid::UUIDgenerate()

    # save EGG object in .RDS format
    saveRDS(EGG, fnEGG)
    if (!silent) {
        print(sprintf("Done: Wrote EGG object to file %s.", fnEGG))
    }

    # ==== WRITE LOG ===========================================================

    if(writeLog) {

        myTitle <- "diffuseV3"

        # Compile function call record
        myCall <- character()
        myCall[1] <- "diffuseV3("
        # ToDo ... update
        myCall[2] <- sprintf("fnAGG = \"%s\", ", fnAGG)
        myCall[3] <- sprintf("fnEGG = \"%s\", ", fnEGG)
        myCall[4] <- sprintf("param = list(")
        myCall[5] <- sprintf("N = %s, ", as.character(N))
        myCall[6] <- sprintf("q = %s, ", as.character(q))
        myCall[7] <- sprintf("seed = %s, ", as.character(seed))
        myCall[8] <- sprintf("useRcpp = %s), ", as.character(useRcpp))
        myCall[9] <- sprintf("silent = %s, ", as.character(silent))
        myCall[10] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")

        # Record progress information
        myNotes <- character()
        myNotes <- c(myNotes, sprintf("Read AGG from file: %s", fnAGG))
        myNotes <- c(myNotes, sprintf("with UUID: %s", attr(AGG, "UUID")))
        myNotes <- c(myNotes, sprintf("wrote EGG to file: %s", fnEGG))
        myNotes <- c(myNotes, sprintf("with UUID: %s", attr(EGG, "UUID")))

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
                 notes = myNotes)
    }

}


# [END]
