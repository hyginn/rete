
consense <- function(sN,minEdgeCC=length(sN),noLog = FALSE,silent = FALSE) {

    #'Find consensus subnetworks
    #'
    #   sN: List of subnetworks resulting after edge removal and
    #'@param sN List of subnetworks
    #'@param minEdgeCC minimum number of subnetworks an edge appears in to be kept
    #'@param noLog if FALSE, save process to an .Rdata log file
    #'@param silent if FALSE, print process to console in real time
    #'@return A list of consensus subnetworks as igraph objects

    #INPUTS:
    #   sN: List of subnetworks
    #   minEdgeCC: minimum number of subnetworks an edge appears in to be kept#   method: method of edge removal
    #   noLog: if FALSE, save process to an .Rdata log file#   delta: threshold for edge removal for given method
    #   silent: if FALSE, print process to console in real time#   EGG: 'Heat Diffusion' annotated AGG, iGraph object.

    #OUTPUT:
    #   cSN: List of subnetworks resulting after edge removal and
    #       subnetwork construction


    #minEdgeCC = 0 #length(outputGraphs)
    numGraphs <- length(outputGraphs)
    sNEdges <- vector(mode="list", length=numGraphs)
    sNVerts <- vector(mode="list", length=numGraphs)
    edgePairsList <- c() #vector(mode="list", length=numGraphs)
    for (i in 1:numGraphs){
        print(i)
        sNEdges <- igraph::as_data_frame(outputGraphs[[i]], what = "edges")
        sNVerts <- igraph::as_data_frame(outputGraphs[[i]], what = "vertices")
        # Find edgepairs
        numEdges <- nrow(sNEdges)
        print(sNEdges)
        edgePairs <- character(length=numEdges) # * 2)
        for (j in 1:numEdges){
            print(edgePairs)
            edgePairs[j] <- paste(sNEdges[j, 'from'], sNEdges[j, 'to'], sep='-')
        }
        edgePairsList <- c(edgePairsList, edgePairs)
    }

    # Weigh edges based on how many subnetworks they occur in
    boolIndex <- as.vector(table(edgePairsList) >= minEdgeCC)
    edgesToKeep <- edgePairsList[boolIndex]

# Keep only edges with weight > minEdgeCC

# Extract connected components

   return(cSN)
}
