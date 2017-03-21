# A value for δ is chosen such that random networks will not result in large subnetworks.

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
#'     Bidirectional edges must only be swapped with each other, therefore the
#'     swapping operation will first be applied to bidirectional edges,
#'     Q * |E_b| times, where E_b is a bidirectional edge. Since all 4 incident
#'     nodes in 2 bidirectional edges are receiving endpoints, one node will be
#'     randomly selected from each edge, and a new edge is formed between the
#'     selected as well as the remaining vertices.
#'     Then Q * |E_u| unidirectional edge swaps take place.
#'     At each unidirectional edge swap, only the receiving endpoints are
#'     swapped e.g. (A → B, C → D) to (A → D, C → B).
#'     In both cases, a swap should NOT result in:
#'     Multiple edges: Possibly result of a bidirectional edge swap.
#'     Self-edges: e.g. (A ↔ B, C ↔ A) to (A ↔ A, B ↔ C)
#'     Disconnected components: Considering (A → B, C → D), if A and B are
#'     in the same component, they should still be in the same component when
#'     (A → D, C → B) This guarantees that the graph is still connected if
#'     the starting graph was also connected.
#'     ### What about (A → B, A → D) to (A → D, A → B)? Will this be
#'     cosidered a swap? (counted towards Q * |E|) ###
#'
#' @section Maintaining Connected Components:
#'     In order to ensure nodes maintain their components during permutation,
#'     the components of the graph is measured initially, and the component
#'     each node belongs to is stored in a separate data structure (hash table?).
#'     After each swap, the new component of the sending endpoint (A in A → B)
#'     is the same as the component the receiving endpoint belongs to. Thus
#'     the component maintenance is acieved in 2 steps (one query and one
#'     comparison) for each receiving node, and thus 4 * Q * |E| steps are
#'     required.
#'
#'     ## Why not use a binary search style approach to find the minimum delta?
#'     ## Can calculation of delta be parallelized?
#' @section Incrementing δ and Lmax: ## How much should delta be incremented for efficieny as well as accuracy?##
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
#' @author Farzan Taj farzantaj@outlook.com
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
    # NOTE: Most comments, only as a confirmation/question


    # Check arguments for validity, generate error if necessary and stop
    # Read EGG
    # Use the helper function ".permute" to create permutations,
    # which takes a copy of EGG, and N, and returns N random
    # graphs as a list (?) of objects


    # If variables are assigned in the local environmet, the helper functions
    # can use them.
    # Generate the list of edges, vertices and vertex count available in EGG
    allEdges <- igraph::get.edgelist(EGG)
    allVertices <- igraph::V(EGG)
    numVertices <- igraph::vcount(EGG)

    # Call permute on EGG, store returned list of graphs as an object
    graphList <- .permute(EGG, N)
    masterDeltaList <- vector("list", length = N)

    # For each graph:

    # Create a vector storing min δ ord for each maxComponentSize
    # (= is used to only store objects in current loop)
    graphDeltas = vector("numeric", Lmax - 1)
    currentDelta = 0.01
    currentGraph = graphList[i]
    currentEdges = igraph::E(currentGraph)

    # Set maxComponentSize to 2 initially, and increase to a max of Lmax.
    maxComponentSize = 2

    # For each Lmax, set δ to 0.01 (?), remove all edges that have weights that
    # are below such δ.
    currentGraph <- igraph::delete.edges(currentGraph, which(currentEdges$weight <= currentDelta))

    # Find the largest 'strongly' connected component:
    largestComponent <- max(igraph::components(currentGraph, mode = "strong")$csize)

    # If largest component has size less than Lmax,
    # continue incrementing δ, else, stop and record δ in the vector, and
    # increase maxComponentSize by 1 (to a max of Lmax)
    if (largestComponent <= maxComponentSize) {
        currentDelta <- currentDelta + 0.01
    } else {
        graphDeltas[i] <- currentDelta
        maxComponentSize = maxComponentSize + 1
        if (maxComponentSize > Lmax) {
            # Must find another way to add to list if parallelization is chosen
            masterDeltaList[[i]] <- graphDeltas
            break
        }
    }
    # (What if the first δ leads to a size equal to maxComponentSize? Should δ
    # be decreased instead? By how much should δ be incremented at each step?)

    # Find the minimum δ for each graph at each maxComponentSize

    # After generating N vectors of size (Lmax - 1), find the median δ for each
    # δ ord (maxComponentSize). i.e. Lmax - 1 median values.
    deltaOrd <- vector("numeric", Lmax - 1)
    for (i in 1:(Lmax - 1)) {
        deltaOrd[i] <- median(sapply(masterDeltaList, FUN = "[", i))
    }

    # Attach the vector as graph metadata to the original EGG and return.
    # (Is this the correct way?)
    attributes(EGG) <- deltaOrd

    # (Can this be parallelized? Since the processes are not interdependent)

    # Write to Console:
    # Update the fraction of N graphs that are done with permutations, (can
    # also update the rough percentage of Q * |E| swaps that are done, but can
    # slow down the process). This is done in .permute.
    # Update the fraction of graphs that are done with finding δ vector (can
    # also indicate the percentatge of maxComponentSize values that are
    # finished per graph, but can again slow the process)

    # Write to log:
    # Log of passed arguments (function call record), processed information
    # (number of edges and vertices, bidirectional and unidirectional)
    # Output object (δ ord vector as metadata), start and finish times
    # (per graph, or for the whole process?)
}



#' Helper function to create N permutations of the EGG graph.
#'
#' @param EGG The given EGG graph from parent function.
#' @param N The number of permuted graphs that must be generated.
#' @param Q The number of edge-swapping permutations to produce one of the
#' permuted networks is Q * |E|.
#'
#' @return List of N igraph objects, permutations of the given EGG graph
.permute <- function(EGG, N, Q) {
    # Calculate the component each vertex belongs to
    ### By subsetting components (somehow)
    initialComponents <- igraph::components(EGG)
    numComponents = initialComponents$no


    # Create a hash table to track components vertecies
    # to. (can a package be used? fastmatch? or just hashed env?)
    # Assign cluster ID to each node
    # Vertex names must be unique
    componentTracker <- new.env(hash = TRUE, size = numVertices)
    for (i in 1:numVertices) {
        # Should $name be used?
       componentTracker[[allVertices$name[i]]] <- initialComponents$membership[i]
    }

    # Get a list of all edges, separate find bidirectional ones
    # Courtesy of Sacha Epskamp,
    E <- t(apply(X = allEdges, MARGIN = 1, FUN = sort))
    bidirectionalEdges <- allEdges[duplicated(E) | duplicated(E, fromLast = TRUE)]
    numBidirectional <- (length(bidirectionalEdges)) / 2
    # Since they are all ordered, we can pick only the first half, and group
    # each pair

    # Better to use a list of tuples to keep track of 'bidirectional' edge pairs
    # as they are simply a pair of unidirectional edges
    bidirPairList <- vector("list", length = (numBidirectional / 2))
    for (i in seq(1, numBidirectional, 2)) {
        bidirPairList[[(i %/% 2) + 1]] <- c(bidirectionalEdges[i], bidirectionalEdges[i + 1])
    }

    ### For each of the N permutated graphs:

    # Create a new copy of the original EGG graph
    currentGraph <- EGG
    ## Swap bidirectional edges first:

    # Randomly pick one from the two nodes involved in a random bidirectional edge,
    # and swap with another such node, ensure to update both edges involved in a
    # bidirectional edge
    i = 0
    # Swap Q * |E_b| times
    while (i < (Q * numBidirectional)) {
        pair <- sample(bidirPairList, size = 2, replace = FALSE)
        # And then change edge participants
        # Must delete mutated edges, and add to graph as new edges (even though
        # the sending nodes are the same)
        # Use add_edges() and delete_edges()
        # Can use igraph::any_multiple() to check if multiple edges are created
        # after each swap, do not increment i if this is the case.

        # Check the new components of all nodes against their previous ones,
        # using the hash structure created above, discard change if components
        # are changed (Does this address the disconnected components issue?)
    }

    # Bidirectional edge swap must not result in multiple edges, thus the new
    # endpoint nodes are checked for the existence of an already existing edge
    # from the newly picked node pairs. Neighbor count will remain the same
    # after a bidirectional edge swap, therefore no need to check.

    # EGG <- random.graph.game(n = 10, p.or.m = 8, type = "gnm", directed = TRUE)


    ## Then swap unidirectional edges:

    # At every unidirectional edge swap, receiving endpoints will have a different
    # number of neighbors, while sending endpoint's node degrees remain the same
    # Therefore the receiving nodes must have another neighbor that targets them
    # also do an edge swap, connecting them to another node.
    # Edge endpoints can be retrieved using:
    # igraph::ends(EGG, es = E(EGG)[vertices])[,2]
    # Therefore edge swaps must be done in two successive pairs:
    # The first pair is selected randomly while the second pair is always the
    # receiving endpoints of the first pair (!!!Doesn't this result in a chain)?

    # igraph::each_edge() can be used to rewire endpoints, can check to see if
    # randomly selected edge is a bidirectional edge (which would be already
    # processed)

    # igraph::keeping_degseq() allows rewiring while preserving "total" degree
    # distribution (of the graph), not sure if it can be used here
    # a <- rewire(graph = EGG, with = keeping_degseq())

    # Function degree() with modes: "in", "out" and "all" can be used to check
    # how many edges point to and from a node.


    ### After both kinds of edge swap (or preparation), each node will be
    # checked to see if their components will be changed. If yes, the swapping
    # process will be aborted, and the operation will not count towards the
    # Q * |E| swaps. (Are none of the nodes allowed to change components?
    # is there a more efficient way of knowing if there are disconnected
    # components are created?)

    # After all the swaps have taken place successfully, the permuted graph
    # will be saved. N permuted graphs are created, and stored in a list.
    # (Can also be saved to drive for future use, as this is computationaly
    # expensive to redo)

    # return list of permuted graphs to caller
}
