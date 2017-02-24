# gGtools.R
#
# Utility functions for gG objects

.df2gG <- function(inFile, arguments, isDirected = TRUE) {
    # Purpose:
    #     Creates a gG object from a dataframe
    #
    # Parameters:
    #     inDF: fileName of the input file
    #     arguments: string - values of calling arguments
    #     isDirected: default TRUE - whether the data frame contains
    #                 directed edges.
    # Details:
    #     The data frame netDF is expected to exist in the calling environment
    #     (i.e. the parent.frame()). igraph expects columns 1 and 2 to contain
    #     vertex names, the remaining columns to contain edge attributes. Edge
    #     attributes will get the name of the column they come from. If
    #     isDirected is true, absent edges are implied to have weight 0 but they
    #     are not explicitly added. If is Directed is false,
    #     igraph::as.directed() expands the graph to have directed edges.
    #     Metadata is attached as a graph attribute.
    # Value:
    #     gG: igraph graph object
    # ToDo:
    #     Check whether existing reverse edges are duplicated by iGraph if
    #     igraph::as.directed() is called on a network.


    # setup metadata
    meta <- list(gGversion = "1.0",
                 logFile = getOption("rete.logFile"),
                 inFile = inFile,
                 args = arguments,
                 date = Sys.Date())

    # create iGraph object
    if (isDirected) {
        gG <- igraph::graph_from_data_frame(get("netDF", parent.frame()),
                                            directed = TRUE)
    } else {
        gG <- igraph::graph_from_data_frame(get("netDF", parent.frame()),
                                            directed = FALSE)
        gG <- igraph::as.directed(gG, mode = "mutual")
    }

    for (name in names(meta)) {
        gG <- igraph::set_graph_attr(gG, name, meta[[name]])
    }

    return(gG)
}

# [END]
