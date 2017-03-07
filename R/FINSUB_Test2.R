#Test 1

#================Test2===============================
folder<-"C:/Users/HPDM1/Documents/CanadaUofT/4thYear/BCB420/ekplektoR/R/"
fPath<-paste(folder,"FINDSUB_Test_Data.R",sep="")
source(fPath) #executing script generating test data. See script for variables generated and how

#delta set to 4, minOrd to 4

igraph::graph_attr(EGG,"delta")<-4 #setting delta parameter (score threshold for edge inclusion)

minOrd <- 4

#We should only get 1 subnetwork this time
outputGraphs<-findsub(method="Leis",EGG,minOrd,noLog = FALSE,silent=FALSE)




#method: method of edge removal; EGG: iGraph object containing
#network with assigned 'influence' to each edge; noLog: do we save a log of the process? ; silent
#if false, print what appears in the process log while the process is in action, as well as additional
#material (e.g. list edges removed, vertices removed)

#Note:

#With a threshold influence of 3.5, LMN5 OPQ6 RST7 should form a network,
#ABC1,DEF2,GHI4,JKL3 should form a network
#When this is lowered to 3, LMN5,OPQ6,RST7,UVW8,and XYZ9 should form a network
#ABC1 DEF2 JKL3 GHI4 from same subnetwork
#Lower further to 2, and the 'discovered subnetwork' will include all vertices, but not all edges
#Lower influence threshold (delta) to 1, and every edge included

#below 2 lines manually constructing expected networks
expectNet1Edges<-BigNetEdges[1:10,]
expectNet1Edges<-expectNet1Edges[order(expectNet1Edges$edgeID),]
expectNet1Verts<-BigNetVertices[1:4,]

#below 3 lines loading iGraph object contents from output into data frames
netwk1<-outputGraphs[[1]]
netwk1Edges<-as_data_frame(netwk1,what="edges")
netwk1Verts<-as_data_frame(netwk1,what="vertices")


#edgeCompar: vector of falses with length(ncol(expectNetXEdges)). If all column in the expected
#dataset are identical to that of the test dataset (netwkXEdges), edgeCompar becomes a vector of
#TRUES with the same length. If this condition is fulfilled, then we know the output for
#the given network is identical to the expected edges.

edgeCompar1<-vector(length=ncol(expectNet1Edges))
if (ncol(expectNet1Edges)==ncol(netwk1Edges) && nrow(expectNet1Edges)==nrow(netwk1Edges)) {
    for (i in 1:ncol(expectNet1Edges)) {
        truRows<-expectNet1Edges[,i]==netwk1Edges[,i]
        sumTruRows<-sum(as.numeric(truRows))
        if (sumTruRows==nrow(expectNet1Edges)) {
            edgeCompar1[i]<-TRUE
        }
    }
}



#vertCompar: because of how unique() sorts things, you can get your vertices out of original order
#Therefore, the aim here is to see that for every vertex/gene score pair in the expected vertex dataframe
#is present in the output vertex data frame (by vertex dataframe I mean dataframe for vertices of
#a given network, expected or produced by the FINDSUB function)
#If that condition is fulfilled, AND the nrow of the expected vertex set is the same as that of the
#output vertex set, then you know that the two data frames of vertices have identical contents, and
# by extension the network's iGraph object's vertices/attributes were added correctly

vertCompar1<-vector(length=nrow(expectNet1Verts))
for (i in 1:nrow(expectNet1Verts)) {
    ind1<-grepl(expectNet1Verts$HGNC_symbol[i],netwk1Verts$name)
    if (sum(as.numeric(ind1))>0) {
        if (netwk1Verts$Gene_Score[i]==expectNet1Verts$Gene_Score[ind1]) {
            vertCompar1[i]<-TRUE
        }
    }
}

#checking the following conditions:
#number of vertices in expected and produced (netwk1) graphs consistent
#every vertex in expected graph is present in produced graph
#every row and column of expected graph reproduced in produced graph
#only one subnetwork produced given input parameters and delta value

context("Check that you have 1 and only one correctly constructed subnetwork")
test_that("1 subnetwork created with expected edges and vertices", {
    expect_equal(nrow(expectNet1Verts),nrow(netwk1Verts))
    expect_equal(sum(as.numeric(vertCompar1)),nrow(expectNet1Verts))
    expect_equal(sum(as.numeric(edgeCompar1)),ncol(expectNet1Edges))
    expect_equal(length(outputGraphs),1)
})

context("Check Heat score")
test_that("Heat scores are as expected", {
    expect_equal(graph_attr(netwk1,"aggHeatScore"),22)
})
