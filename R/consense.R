HGNCsymb<-c("ABC1","DEF2","GHI4","JKL3","LMN5","OPQ6","RST7","UVW8","XYZ9")
VertScores<-c(10,2,4,6,7,8,2,5,9)

BigNetVertices<-data.frame(HGNC_symbol=HGNCsymb,Gene_Score=VertScores)
From<-c("ABC1","ABC1","DEF2","JKL3","DEF2","DEF2","JKL3","GHI4","JKL3","GHI4","GHI4","GHI4","LMN5",
        "XYZ9","LMN5","OPQ6","OPQ6","UVW8","OPQ6","RST7","UVW8","XYZ9")
To<-c("DEF2","JKL3","ABC1","ABC1","JKL3","GHI4","DEF2","DEF2","GHI4","JKL3","LMN5","XYZ9","GHI4",
      "GHI4","OPQ6","LMN5","UVW8","OPQ6","RST7","OPQ6","XYZ9","UVW8")
edgeID<-character(length=length(From))

for (i in 1:length(edgeID)) {
    edgeID[i]<-paste("BNE_",as.character(i),sep="")
}



#Influence vector created to assign influence levels to directed edges
Influence<-numeric(length(From))
Influence[1:length(Influence)]=5
Influence[(From=="GHI4")&(To=="LMN5")]=1
Influence[(From=="LMN5")&(To=="GHI4")]=2
Influence[(From=="GHI4")&(To=="XYZ9")]=1.5
Influence[(From=="XYZ9")&(To=="GHI4")]=2.5
Influence[(From=="OPQ6")&(To=="UVW8")]=2
Influence[(From=="UVW8")&(To=="OPQ6")]=3
Influence[(From=="UVW8")&(To=="XYZ9")]=3.5
Influence[(From=="XYZ9")&(To=="UVW8")]=3.5
#With a threshold influence of 3.5, LMN5 OPQ6 RST7 should form a network,
#ABC1,DEF2,GHI4,JKL3 should form a network
#When this is lowered to 3, LMN5,OPQ6,RST7,UVW8,and XYZ9 should form a network
#ABC1 DEF2 JKL3 GHI4 from same subnetwork
#Lower further to 2, and the 'discovered subnetwork' will include all vertices, but not all edges
#Lower influence threshold (delta) to 1, and every edge included

BigNetEdges<-data.frame(from=From,to=To,Influence=Influence,edgeID=edgeID) #edge data frame created for EGG
EGG<-igraph::graph_from_data_frame(BigNetEdges,directed=TRUE,vertices=BigNetVertices) #setup EGG graph
igraph::graph_attr(EGG, "EGGversion") <- "1.0"  #set EGG version attribute
#Here lies the border of the common variable setup.


igraph::graph_attr(EGG,"delta")<-4 #setting delta parameter (score threshold for edge inclusion)

minOrd <- 2  #setting minOrd to 2

#below, applying findSub function
outputGraphs<-findSub(method="Leis",EGG,minOrd,noLog=FALSE,silent=FALSE)
#------------------------------------------------------------------------#


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
