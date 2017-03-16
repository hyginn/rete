# gGtools.R
#
# Utility functions for gG objects

.df2gG <- function(netDF, isDirected = TRUE, simplify = TRUE) {
    # Purpose:
    #     Creates a gG object from the dataframe netDF
    #
    # Parameters:
    #     netDF: a dataframe with two columns of vertex IDs and one column
    #            of edge weights.
    #     isDirected: default TRUE - whether the data frame contains
    #                 directed edges.
    #     simplify: default TRUE - whether or not to collapse multiple edges
    # Details:
    #     igraph expects columns 1 and 2 of netDF to contain vertex names, the
    #     remaining columns to contain edge attributes. Edge attributes will get
    #     the name of the column they come from. If isDirected is true, absent
    #     edges are implied to have weight 0 but they are not explicitly added.
    #     If is Directed is false, igraph::as.directed() expands the graph to
    #     have directed edges. If simplify is true multiple edges and loops are
    #     collapsed and the max() of the weights is the attribute of the
    #     combined edge.
    #
    # gG objects are defined to be weighted, directed, simple graphs
    # the parameters isDirected and simplify are provided merely for
    # development purposes and should not be changed.

    # Value:
    #     gG: igraph graph object with metadata attached as attributes

    # ==== SETUP METADATA ======================================================


    meta <- list(type = "gG",
                 version = "1.0",
                 UUID = uuid::UUIDgenerate())

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
    # ToDo - report if edges were simplified away since this
    # may give us less than the requested number xN of edges.


    # ==== ATTACH METADATA =====================================================
    for (name in names(meta)) {
        attr(gG, name) <- meta[[name]]
    }

    return(gG)
}

# [END]
