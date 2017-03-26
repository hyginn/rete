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
#'     It is sufficient to find all articulation points and identify edges that
#'     connect to them (incoming and outgoing), and during a swap, an edge
#'     between two adjacent articulation points is not to be swapped with an
#'     incoming edge to another articulation point. Therefore EGG's articulation
#'     points are initially discovered by Tarjan's algorithm and
#'     To ensure nodes maintain their components during permutation,
#'     the size and distribution of the components of the graph is measured
#'     initially, and and only
#'     the largest connected component is permutated.
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
#'     the graph for each Lmax starting from 2 to the Lmax given as an argument.
#'     On each such run, and a small value for δ is chosen (0.01) , and
#'     the strongly connected component discovery is done. If the size of the
#'     biggest component after such run is smaller than Lmax, increment δ by
#'     (0.01), and redo component disovery. As soon as the largest component has
#'     size Lmax (for current run), δ is recorded, and Lmax is incremented by 1.
#'     Once current Lmax is the same as the given Lmax, run the procedure on the
#'     next graph.
#'     The above procedure is done on each of the N graphs, Lmax - 1 times.
#'
#'
#'
#'
#' @param EGG An EGG object - the input heat-equilibrated gene score graph.
#' @param Lmax	Maximum order for largest connected component in a random
#' network. See details.
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
                   N = 100,
                   Q = 100,
                   silent = FALSE,
                   writeLog = TRUE) {
    # ==== VALIDATIONS =========================================================

    # Check if EGG is given
    if (missing(EGG)) {
        stop(paste0("Missing input EGG, with no defaults."))
    }
    # General parameter checks
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

    # Create a list of deltas for each permuted graph
    masterDeltaList <- vector("list", length = N)

    for (i in 1:N) {
        # Permute graph with given Q
        cat(sprintf("Creating permutation %d of %d\n", i, N))
        currentGraph <- .permuteGraph(EGG, Q)
        cat(sprintf("\nCalculating deltas for graph %d of %d\n", i, N))

        # Set vector to hold deltas, initial delta, maxComponentSize (ord)
        graphDeltas <- vector("numeric", Lmax - 1)
        currentDelta <-  0.0005
        maxComponentSize <-  2
        j <- 1
        k <- 0
        testGraph <- currentGraph

        # For each value of maxComponentSize (ord):
        while (maxComponentSize <= Lmax) {

            # Find edges that have weights lower than current delta and delete
            edgesToRemove = igraph::E(testGraph)[E(testGraph)$weight < currentDelta]
            testGraph = igraph::delete.edges(testGraph, edgesToRemove)

            # Find the size of the largest strongly connected component
            currentLargest <-
                max(igraph::components(testGraph, "strong")$csize)

            # If largest component has size less than Lmax, continue incrementing δ
            # else, stop and record δ in the vector, and increase maxComponentSize
            # by 1
            if (currentLargest >= maxComponentSize) {
                currentDelta <- currentDelta + 0.0005
            } else {
                .pBar(k, Lmax - 1, 80)
                # Reset and increase ord
                graphDeltas[j] <- currentDelta
                maxComponentSize = maxComponentSize + 1
                currentDelta <- 0.0005
                j <- j + 1
                k <- k + 1
                testGraph <- currentGraph
            }
        }

        # Append this graph's deltas to the master list of deltas
        masterDeltaList[[i]] <- graphDeltas
    }

    # By how much should δ be incremented at each step?)

    # After generating N vectors of size (Lmax - 1), find the median δ for each
    # δ ord (maxComponentSize). i.e. Lmax - 1 median values.
    deltaOrd <- vector("numeric", Lmax - 1)
    print("Calculating median deltas")
    for (i in 1:(Lmax - 1)) {
        .pBar(i, Lmax - 1, Lmax - 1)
        deltaOrd[i] <- median(sapply(masterDeltaList, FUN = "[", i))
    }

    # Attach the vector as graph metadata to the original EGG and return.

    igraph::graph_attr(EGG, "delta") <- deltaOrd
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
    # ==== Prepare largest component for permutation ===========================

    # Get all components in the EGG
    # allComponents <- igraph::components(EGG)

    # Select the subgraph containing the largest component
    mainComponent <- decompose(EGG, "weak", 1)[[1]]

    # tempTest <- mainComponent

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
        assign(x = as.character(bidirectionalEdges[i]), value = bidirectionalEdges[i + (nBi / 2)], envir = biDirs)
        assign(x = as.character(bidirectionalEdges[i + (nBi / 2)]), value = bidirectionalEdges[i], envir = biDirs)

        # biDirs[[as.character(bidirectionalEdges[i])]] <- bidirectionalEdges[i + (nBi / 2)]
        # biDirs[[as.character(bidirectionalEdges[i + (nBi / 2)])]] <- bidirectionalEdges[i]
    }

    # ==== Swap edges ==========================================================
    # NOTE: Edge IDs must be hardcoded as attributes, as the actual IDs change
    # after each swap and cannot be tracked.
    tempEdgeVs =  unlist(lapply(1:nEdges, function(i)
        allEdges[i, ]))
    tempEdgeIDs = igraph::get.edge.ids(mainComponent, vp = tempEdgeVs, directed = TRUE)
    mainComponent <-
        igraph::set.edge.attribute(mainComponent, name = "ID", value = tempEdgeIDs)

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
                    curVP1 = as.character(pair[1, ])
                    curVP2 = as.character(pair[2, ])
                    N1 = unlist(lapply(curVP1, function(x)
                        igraph::neighbors(mainComponent, v = x, mode = "all")))
                    N2 = unlist(lapply(curVP2, function(x)
                        igraph::neighbors(mainComponent, v = x, mode = "all")))

                    if (any(N1 %in% N2)) {
                        # Skip iteration, as this would result in a multi-edge
                        next
                    }
                    # ==== Bidirectional edge swap =================================================

                    # Get edge weights from subgraph, prepare as attributes
                    # TODO: There may be a bug here...
                    eIDs = igraph::E(graph = mainComponent,
                                     P = c(pair[1, ], pair[2, ], rev(pair[1, ]), rev(pair[2, ])))$ID
                    eWeights = igraph::E(mainComponent)$weight[E(mainComponent)$ID %in% eIDs]
                    attribs = list(weight = eWeights, ID = eIDs)

                    # Sample a vertex from incident vertices, pick remainders
                    v1 = sapply(lapply(list(pair[1, ], pair[2, ]), sample), function(x)
                        x[1])
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
                        igraph::delete.edges(mainComponent, c(del1, del2, del3, del4))
                    mainComponent <-
                        igraph::add.edges(mainComponent,
                                          c(v1, v2, revV1, revV2),
                                          attr = attribs)

                    # Update biDirs
                    biDirs[[as.character(v1[1])]] <- v1[2]
                    biDirs[[as.character(v1[2])]] <- v1[1]

                    biDirs[[as.character(v2[1])]] <- v2[2]
                    biDirs[[as.character(v2[2])]] <- v2[1]

                    # Update i
                    i <- i + 1
                }
            } else {
                # ==== Unidirectional edge swap ================================================

                # Get edge weights from subgraph, prepare as attributes
                eIDs = igraph::E(mainComponent, c(pair[1, ], pair[2, ]))$ID
                eWeights = igraph::E(mainComponent)$weight[E(mainComponent)$ID %in% eIDs]

                attribs = list(weight = eWeights, ID = eIDs)

                # Swap endpoints
                v1 = c(pair[1, 1], pair[2, 2])
                v2 = c(pair[2, 1], pair[1, 2])

                # Form new edges, delete old edges
                del1 = paste0(pair[1, ], collapse = "|")
                del2 = paste0(pair[2, ], collapse = "|")
                mainComponent <-
                    igraph::delete.edges(mainComponent, c(del1, del2))
                mainComponent <-
                    igraph::add.edges(mainComponent, c(v1, v2), attr = attribs)

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

# [END]
