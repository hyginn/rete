# gGtools.R
#
# Utility functions for gG objects

.df2gG <- function(inFile, call, isDirected = TRUE, simplify = TRUE) {
    # Purpose:
    #     Creates a gG object from a dataframe
    #
    # Parameters:
    #     inDF: fileName of the input file
    #     arguments: string - values of calling arguments
    #     isDirected: default TRUE - whether the data frame contains
    #                 directed edges.
    #     simplify: default TRUE - whether or not to collapse multiple edges
    # Details:
    #     The data frame netDF is expected to exist in the calling environment
    #     (i.e. the parent.frame()). igraph expects columns 1 and 2 to contain
    #     vertex names, the remaining columns to contain edge attributes. Edge
    #     attributes will get the name of the column they come from. If
    #     isDirected is true, absent edges are implied to have weight 0 but they
    #     are not explicitly added. If is Directed is false,
    #     igraph::as.directed() expands the graph to have directed edges.
    #     Metadata is attached as object attributes. If simplify is true
    #     multiple edges and loops are collapsed and the max() of the weights
    #     is the attribute of the combined edge.
    # Value:
    #     gG: igraph graph object with metadata attached as object attributes
    # ToDo:
    #     Check whether existing reverse edges are duplicated by iGraph if
    #     igraph::as.directed() is called on a network.


    # ==== SETUP METADATA ======================================================
    meta <- list(version = "gG 1.0",
                 UUID = "12345",
                 inFile = inFile,
                 call = call,
                 time = Sys.time())

    # ==== CREATE IGRAPH GRAPH =================================================
    if (isDirected) {
        gG <- igraph::graph_from_data_frame(get("netDF", parent.frame()),
                                            directed = TRUE)
    } else {
        gG <- igraph::graph_from_data_frame(get("netDF", parent.frame()),
                                            directed = FALSE)
        gG <- igraph::as.directed(gG, mode = "mutual")
    }

    # ==== SIMPLIFY GRAPH ======================================================
    if (simplify) {
        gG <- igraph::simplify(gG,
                               remove.multiple = TRUE,
                               remove.loops = TRUE,
                               edge.attr.comb = "max")
    }
    # ToDo - post log message if edges were simplified away since this
    # may give us less than the requested number xN of edges.


    # ==== ATTACH METADATA =====================================================
    for (name in names(meta)) {
        if (! is.null(attr(gG, name))) {
            stop(sprintf("Error: can't overwrite existing attribute %s.", name))
        }
        attr(gG, name) <- meta[[name]]
    }

    return(gG)
}

# [END]
