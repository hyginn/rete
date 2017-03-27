# edgeThresh.R

#' Define edge-removal threshold parameter δ for each target subnet-order, from
#' random permutations of the input EGG object
#'
#' \code{THRESH} creates N permutations of the input graph EGG by randomly
#' choosing 2 edges and swapping their target endpoints (i.e. B in the
#' edge A → B). Finds the value δ (edge weight) for each random graph such that
#' the maximum size of any component in the graph is Lmax, and returns the
#' median of N δ values.
#'
#' @section Swapping Edges:
#'     Bidirectional edges must only be swapped with each other, Since all 4
#'     incident nodes in 2 bidirectional edges are receiving endpoints, one node
#'     will be randomly selected from each edge, and a new edge is formed
#'     between the selected as well as the remaining vertices.
#'     At each unidirectional edge swap, only the receiving endpoints are
#'     swapped e.g. (A → B, C → D) to (A → D, C → B).
#'     In both cases, a swap should NOT result in:
#'
#'     "Multiple edges": Possibly result of multiple bidirectional edge swaps:
#'             A                             A
#'          //   \\              →                 //
#'         B   →   C                     B  ↔ →  C
#'
#'     OR:
#'     Where an edge belonging to vertex (A) which already has an outgoing edge to
#'     another vertex (B) is switched with another incoming edge to vertex B,
#'     resulting in (A) having two outgoing edges to (B).
#'
#'        ( )     ( )                  ( )  ←  ( )
#'         ↑       ↓             →
#'        (A)  →  (B)                  (A)   ⇉  (B)
#'
#'     "Self-edges": e.g. (A ↔ B, C ↔ A) to (A ↔ A, B ↔ C)
#'     Where in the case of:
#'     'Unidirectional edges': Precedent of one edge is the antecedent of the
#'     other edge to be swapped.
#'     'Bidirectional edges': The edge pair to be swapped involves less than 4
#'     distinct nodes.
#'     "Disconnected components": Considering (A → B, C → D), since A and B are
#'     in the same component, they should still be in the same component when
#'     (A → D, C → B). This guarantees that the component remains connected.
#'
#'     Note that (A → B, A → D) to (A → D, A → B) is not considered a swap, as
#'     an edge swap has to involve 4 distinct vertices. This will also avoid
#'     multiple edges and self-edges, as they involve less than 4 distinct
#'     vertices.
#'
#' @section Maintaining Connected Components:
#'     Note: In order to minimize computation, and since the target outcome is
#'     to measure order of the largest connected component, only the largest
#'     connected component of the input is then used for the computation in case
#'     the input has more than one component.
#'
#'     To ensure nodes maintain their components during permutation,
#'     the size and distribution of the components of the graph is measured
#'     initially, and and only the largest connected component is permutated.
#'
#'     Then the strongly connected components within the largest component are
#'     found.
#'     After each swap, the new component of the sending endpoint (A in A → B)
#'     is the same as the component the receiving endpoint belongs to. Thus
#'     the component maintenance is achieved in 2 steps (one query and one
#'     comparison) for each receiving node, and thus 4 * Q * |E| steps are
#'     required.
#'
#' @section Incrementing δ and Lmax:
#'     After creating the permutated graph, the following procedure is run on
#'     the graph for each Lmax starting from Lmin to the Lmax arguments.
#'     Initially, δ values are chosen by finding min(|E|, 100) quantiles of
#'     edge weight distribution, and used to delete edges from the permutation.
#'     As soon as the largest component has size (ord) between Lmin and Lmax,
#'     If the current permutation's ord values is less than Lmin (without
#'     finding δ), the delta is skipped, and interpolated later on.
#'
#' @param EGG An EGG object - the input heat-equilibrated gene score graph.
#' @param Lmax	Maximum order for the largest connected component in a random
#' network. See details.
#' @param Lmin Minimum order for the largest connected component. Defaults to 2.
#' @param N	The number of permuted networks to produce for choosing a δ
#' cutoff.
#' @param Q	The number of edge-swapping permutations to produce one of the
#' permuted networks is Q * |E|.
#' @param silent Boolean. Default: FALSE. Whether or not to write
#' progress information to console.
#' @param writeLog Boolean. Default: TRUE Whether or not to log results.
#'
#' @return The Equilibrated Gene Graph EGG with δ values as additional
#' graph metadata.
#'
#' @references
#' @seealso
#' @examples
#'
#' @export
THRESH <- function(EGG,
                   Lmax = max(getOptions("rete.minSubnetOrder")),
                   Lmin = 2,
                   N = 100,
                   Q = 100,
                   silent = FALSE,
                   writeLog = TRUE) {
    # ==== VALIDATIONS =========================================================

    # Check if EGG is given
    if (missing(EGG)) {
        stop(paste0("Missing input EGG, with no defaults."))
    }
    # General parameter checks (TODO)
    # cR <- character()
    # cR <-
    #     c(cR, .checkArgs(EGG,         like = getOptions("rete.EGGprototype")))
    # cR <-
    #     c(cR, .checkArgs(Lmax,        like = numeric(1), checkSize = TRUE))
    # cR <-
    #     c(cR, .checkArgs(N,           like = numeric(1), checkSize = TRUE))
    # cR <-
    #     c(cR, .checkArgs(Q,           like = numeric(1), checkSize = TRUE))
    # cR <-
    #     c(cR, .checkArgs(silent,      like = logical(1), checkSize = TRUE))
    # cR <-
    #     c(cR, .checkArgs(writeLog,    like = logical(1), checkSize = TRUE))

    # if (length(cR) > 0) {
    #     stop(cR)
    # }
    # Lmax cannot be more than graph size
    # Validate the Lmax parameter
    if (!missing(Lmax) &&
        (Lmax > (numVertices <- igraph::vcount(EGG)))) {
        stop(paste0("Lmax cannot be more than graph size."))
    }
    if (!missing(Lmax) && (Lmax < 2)) {
        stop(paste0("Lmax cannot be less than 2."))
    }
    # Validate that the starting graph's largest connected component (weakly)
    if (max(igraph::components(EGG, "weak")$csize) < Lmax) {
        stop(paste0(
            "Largest connected component in given graph has size
            smaller than Lmax."
        ))
    }
    # ((Report if graph already has delta attributes?))
    # ((Is attribute check enough or are there other attributes?))
    if (!is.null(igraph::graph.attributes(EGG)$delta)) {
        stop(paste0("Graph already has attributes, remove and retry."))
    }
    # ((Should extremely large N and Q be tested for? Maybe a warning should be
    # printed))

    # ==== Permute and find deltas==============================================

    # Use the helper function ".permuteGraph" to create permutations,
    # which takes a copy of EGG, and returns a random
    # graph as an igraph object

    # Create a data frame to hold each permuted graph's deltas
    masterDataFrame <-
        data.frame(matrix(
            vector(),
            nrow = N,
            ncol = (Lmax - Lmin + 1),
            dimnames = list(N = 1:N, ords = (sapply(Lmin:Lmax, function(x)
                paste0("ORD", x))))
        ))

    # Number of steps/increments, and their values
    nQuant <- min(igraph::ecount(EGG), 100)
    edgeWeights <- igraph::E(EGG)$weight
    stepDeltas <- stats::quantile(edgeWeights, seq(0, 1, (1 / nQuant)))

    # Permute and calculate deltas
    for (i in 1:N) {
        # Permute graph with given Q
        cat(sprintf("\nCreating permutation %d of %d\n", i, N))
        currentGraph <- .permuteGraph(EGG, Q)

        cat(sprintf("\nCalculating deltas for graph %d of %d\n", i, N))

        for (j in 1:nQuant) {
            .pBar(j, nQuant, 80)
            # Find edges that have weights lower than current delta and delete
            edgesToRemove = igraph::E(currentGraph)[igraph::E(currentGraph)$weight < stepDeltas[j]]
            currentGraph = igraph::delete.edges(currentGraph, edgesToRemove)

            # Find the size of the largest strongly connected component, ord
            ord <- max(igraph::components(currentGraph, "strong")$csize)

            # If largest component has size less than Lmax, continue incrementing δ
            # else, stop and record δ in the vector, and increase ord
            # by 1
            if (ord < Lmin) {
                break
            } else if (ord < Lmax) {
                # Record ord at correct position
                masterDataFrame[i, ord] <- stepDeltas[j]
            }
        }
    }

    # Interpolate NAs, row wise
    deltaList = lapply(1:(Lmax - Lmin + 1), function(i)
        .interpolateNA(unlist(masterDataFrame[i,]), 0, 1))

    # Get means, "column" wise
    deltaOrds = sapply(1:(Lmax - Lmin + 1), function(i) {
        # Get the i'th element from each sublist of deltaList
        temp = sapply(deltaList, "[[", i)
        # Return means
        return(mean(temp))
    })

    # Attach the vector as graph metadata to the original EGG and return.
    igraph::graph_attr(EGG, "delta") <- deltaOrds
    return(EGG)
}

# TODO Write to log:
# Log of passed arguments (function call record), processed information
# (number of edges and vertices, bidirectional and unidirectional)
# Output object (δ ord vector as metadata), start and finish times
# (per graph, or for the whole process?)




# ==== Graph Permutation =======================================================
#' Helper function to create a valid permutation of the EGG graph.
#' Only permutates the largest connected component
#' @param EGG The given EGG graph from parent function.
#' @param Q The number of edge-swapping permutations to produce the
#' permuted network is Q * |E|.
#'
#' @return Permuted largest component of EGG
.permuteGraph <- function(EGG, Q) {

    # Parameter check
    cR <- character()
    cR <- c(cR, .checkArgs(EGG, like = getOptions("rete.EGGprototype")))
    cR <- c(cR, .checkArgs(Q ,like = numeric(1), checkSize = TRUE))
    if (length(cR) > 0) {
        stop(cR)
    }

    # ==== Prepare largest component for permutation ===========================

    # Select the subgraph containing the largest component
    mainComponent <- igraph::decompose(EGG, "weak", 1)[[1]]

    # Gather component information
    allEdges <-  igraph::get.edgelist(mainComponent)
    allVertices <- igraph::V(mainComponent)
    nVertices <- igraph::vcount(mainComponent)
    nEdges <- igraph::ecount(mainComponent)

    # Recognize bidirectional edges
    # Courtesy of Sacha Epskamp,
    tempEdges = t(apply(
        X = allEdges,
        MARGIN = 1,
        FUN = sort
    ))
    bidirectionalEdges = allEdges[duplicated(tempEdges) |
                                      duplicated(tempEdges, fromLast = TRUE)]
    nBi = length(bidirectionalEdges)

    # Add bidirectional edges to a hashed environment for query (add both
    # vertices)
    # ((Should we use fastmatch instead? Or named lists?))
    biDirs <- new.env(size = nBi)
    for (i in 1:(nBi / 2)) {
        assign(x = as.character(bidirectionalEdges[i]),
               value = bidirectionalEdges[i + (nBi / 2)], envir = biDirs)
        assign(x = as.character(bidirectionalEdges[i + (nBi / 2)]),
               value = bidirectionalEdges[i], envir = biDirs)
    }

    # ==== Swap edges ==========================================================
    # NOTE: Edge IDs must be hardcoded as attributes, as the actual IDs change
    # after each swap and cannot be tracked.
    tempEdgeVs =  unlist(lapply(1:nEdges, function(i)
        allEdges[i, ]))
    tempEdgeIDs = igraph::get.edge.ids(mainComponent, vp = tempEdgeVs,
                                       directed = TRUE)
    mainComponent <- igraph::set.edge.attribute(mainComponent, name = "ID",
                                   value = tempEdgeIDs)

    # Swap Q * |E| times
    i <- 0
    nSwap <- (Q * nEdges)

    while (i < nSwap) {
        # Update progress bar
        .pBar(i, nSwap, 80)

        # Must use current edges in the subgraph:
        allEdges <- igraph::get.edgelist(mainComponent)

        # Randomly select two edges
        pair <- allEdges[sample(1:nEdges, size = 2, replace = FALSE), ]

        # Check if the swap will result in a multiple edge
        # e.g. where for (A,B),(C,D), (A,D) OR (C,B) already exist
        if (any(utils::tail(duplicated(rbind(
            allEdges, c(pair[1, 1], pair[2, 2])
        )))) ||
        any(utils::tail(duplicated(rbind(
            allEdges, c(pair[2, 1], pair[1, 2])
        ))))) {
            next
        }

        # Check if the vertices are all distinct
        if (length(unique(as.vector(pair))) == 4) {
            # Check if any of these edges are bidirectional
            check <-
                mget(as.character(c(pair[1,], pair[2,])), biDirs, ifnotfound = NA)
            if (any(!is.na(check))) {
                if (!all(!is.na(check))) {
                    # Not all are bidirectional, skip (TODO biDirs have a miniscule
                    # chance of swapping this way...)
                    next
                } else {
                    # Check if any of the vertices in the first pair has an edge
                    # to the second pair
                    curVP1 = as.character(pair[1,])
                    curVP2 = as.character(pair[2,])
                    N1 = unlist(lapply(curVP1, function(x)
                        igraph::neighbors(mainComponent, v = x, mode = "all")))
                    N2 = unlist(lapply(curVP2, function(x)
                        igraph::neighbors(mainComponent, v = x, mode = "all")))

                    if (any(N1 %in% N2)) {
                        # Skip iteration, as this would result in a multi-edge
                        next
                    }
                    # ==== Bidirectional edge swap =============================

                    # Get edge weights from subgraph, prepare as attributes
                    # TODO: There may be a bug here...
                    eIDs = igraph::E(graph = mainComponent,
                                     P = c(pair[1,], pair[2,], rev(pair[1,]),
                                           rev(pair[2,])))$ID
                    eWeights = igraph::E(mainComponent)$weight[igraph::E(mainComponent)$ID %in% eIDs]
                    attribs = list(weight = eWeights, ID = eIDs)

                    # Sample a vertex from incident vertices, pick remainders
                    v1 = sapply(lapply(list(pair[1, ], pair[2, ]), sample),
                                function(x) x[1])
                    vVector = as.vector(pair)
                    v2 = vVector[!(vVector %in% v1)]
                    revV1 = rev(v1)
                    revV2 = rev(v2)

                    # Form new edges, delete old edges
                    del1 = paste0(pair[1, ], collapse = "|")
                    del2 = paste0(pair[2, ], collapse = "|")
                    del3 = paste0(rev(pair[1, ]), collapse = "|")
                    del4 = paste0(rev(pair[2, ]), collapse = "|")
                    mainComponent <-
                        igraph::delete.edges(mainComponent, c(del1, del2, del3,
                                                              del4))
                    mainComponent <-
                        igraph::add.edges(mainComponent,
                                          c(v1, v2, revV1, revV2),
                                          attr = attribs)

                    # Update biDirs
                    assign(x = as.character(v1[1]), value = v1[2], biDirs)
                    assign(x = as.character(v1[2]), value = v1[1], biDirs)

                    assign(x = as.character(v2[1]), value = v2[2], biDirs)
                    assign(x = as.character(v2[2]), value = v2[1], biDirs)

                    # Update i
                    i <- i + 1
                }
            } else {
                # ==== Unidirectional edge swap ================================

                # Get edge weights from subgraph, prepare as attributes
                eIDs = igraph::E(mainComponent, c(pair[1, ], pair[2, ]))$ID
                eWeights = igraph::E(mainComponent)$weight[igraph::E(mainComponent)$ID %in% eIDs]

                attribs = list(weight = eWeights, ID = eIDs)

                # Swap endpoints
                v1 = c(pair[1, 1], pair[2, 2])
                v2 = c(pair[2, 1], pair[1, 2])

                # Form new edges, delete old edges
                del1 = paste0(pair[1, ], collapse = "|")
                del2 = paste0(pair[2, ], collapse = "|")
                mainComponent <- igraph::delete.edges(mainComponent,
                                                      c(del1, del2))
                mainComponent <- igraph::add.edges(mainComponent,
                                                   c(v1, v2), attr = attribs)

                # Update i
                i <- i + 1
            }
        } else {
            # Selected vertices are not all distinct, skip
            next
        }
    }

    return(mainComponent)

}

#' interpolates NA values in a vector by interpolateing from two
#' bracketing non-missing values. If is.na(v[1]), the left-bracket is
#' vMin, if is.na(v[length(v)]), vMax is the right bracket.
#' @param v a numeric vector
#' @param vMin numeric
#' @param vMax numeric
#'
#' @return a numeric vector without NAs
.interpolateNA <- function(v, vMin, vMax) {

    if (length(v) == 0) { return(numeric()) }

    # ensure left and rightmost value of v is not NA
    iMin <- 1                                 # index of first element to return
    iMax <- length(v)                         # index of last element to return
    if (is.na(v[1])) {
        v <- c(vMin, v)
        iMin <- iMin + 1
        iMax <- iMax + 1
    }
    if (is.na(v[length(v)])) {
        v <- c(v, vMax)
    }

    # analyze the vector with rle()
    rleV <- rle(is.na(v))                       # run length encoding of is.na()
    rleV$last <- cumsum(rleV$lengths)           # last index of every run
    rleV$first <- rleV$last - rleV$lengths + 1  # first index of every run

    naRun <- which(rleV$values)                 # indices of NA-value runs
    for (i in naRun) {
        iLowerBracket <- rleV$first[i] - 1
        iUpperBracket <- rleV$last[i] + 1
        newV <- seq(v[iLowerBracket],           # value in lower-bracket element
                    v[iUpperBracket],           # value in upper-bracket element
                    length.out = iUpperBracket - iLowerBracket + 1)
        v[iLowerBracket:iUpperBracket] <- newV  # replace values with
                                                # interpolated sequence
    }

    return(v[iMin:iMax])
}

# TEST
# EGG <- igraph::random.graph.game(1000, p.or.m = 3000, type = "gnm", directed = TRUE)
# E(EGG)$weight <- runif(n = igraph::ecount(EGG), min = 0, max = 1)
# (DONE <- THRESH(EGG = EGG, Lmax = 20, N = 1, Q = 1))

# [END]
