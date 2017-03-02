#Test 4

#================Test4===============================
folder<-"C:/Users/HPDM1/Documents/CanadaUofT/4thYear/BCB420/ekplektoR/R/"
fPath<-paste(folder,"FINDSUB_Test_Data.R",sep="")
source(fPath) #executing script generating test data. See script for variables generated, network structure



#Here, we will be assigning an edge influence attribute of 5 to edge GHI4->XYZ9. Therefore,
#ABC1 DEF2 JKL3 GHI4 XYZ9 and UVw8 should form 1 large subnetwork, and LMN5, OPQ6, RST7 should form
#a subnetwork of 3 nodes.

#The main purpose of this test is to show that construction of subnetworks can involve grouping
#edges that are not grouped together i.e. consecutive in the list

Influence[(From=="UVW8")&(To=="XYZ9")]=5
Influence[(From=="XYZ9")&(To=="UVW8")]=5
Influence[(From=="GHI4")&(To=="XYZ9")]=5

BigNetEdges<-data.frame(from=From,to=To,Influence=Influence,edgeID=edgeID) #remaking the data frame
inputGraph<-graph_from_data_frame(BigNetEdges,directed=TRUE,vertices=BigNetVertices) #and inputGraph

outputGraphs<-FINDSUBv1_2(MethType="Leis",Thresh=4,inputGraph,minOrd=3,logResults = TRUE,silent=FALSE)
#MethType: method of edge removal; Thresh: influence threshold; inputgraph: iGraph object containing
#network with assigned 'influence' to each edge; logResults: do we save a log of the process? ; silent
#if false, print what appears in the process log while the process is in action, as well as additional
#material (e.g. list edges removed, vertices removed)



#below 6 lines: set up expected edges and vertices
expectNet1Edges<-BigNetEdges[c(1:10,12,21,22),]
expectNet1Verts<-BigNetVertices[c(1:4,8:9),]
expectNet2Edges<-BigNetEdges[c(15:16,19:20),]
expectNet2Verts<-BigNetVertices[5:7,]
#below 9 lines: set up edges and vertices from output as data frames
netwk1<-outputGraphs[[1]]
netwk1Edges<-as_data_frame(netwk1,what="edges")
netwk1Verts<-as_data_frame(netwk1,what="vertices")
netwk2<-outputGraphs[[2]]
netwk2Edges<-as_data_frame(netwk2,what="edges")
netwk2Verts<-as_data_frame(netwk2,what="vertices")


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

edgeCompar2<-vector(length=ncol(expectNet2Edges))
if (ncol(expectNet2Edges)==ncol(netwk2Edges) && nrow(expectNet2Edges)==nrow(netwk2Edges)) {
    for (i in 1:ncol(expectNet2Edges)) {
        truRows<-expectNet2Edges[,i]==netwk2Edges[,i]
        sumTruRows<-sum(as.numeric(truRows))
        if (sumTruRows==nrow(expectNet2Edges)) {
            edgeCompar2[i]<-TRUE
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
vertCompar2<-vector(length=nrow(expectNet2Verts))
for (i in 1:nrow(expectNet2Verts)) {
    ind2<-grepl(expectNet2Verts$HGNC_symbol[i],netwk2Verts$name)
    if (sum(as.numeric(ind1))>0) {
        if (netwk2Verts$Gene_Score[i]==expectNet2Verts$Gene_Score[ind2]) {
            vertCompar2[i]<-TRUE
        }
    }
}

context("Check that you have 3 subnetworks")
test_that("3 subnetworks created with expected edges and vertices", {
    expect_equal(nrow(expectNet1Verts),nrow(netwk1Verts))
    expect_equal(nrow(expectNet2Verts),nrow(netwk2Verts))
    expect_equal(sum(as.numeric(vertCompar1)),nrow(expectNet1Verts))
    expect_equal(sum(as.numeric(vertCompar2)),nrow(expectNet2Verts))
    expect_equal(sum(as.numeric(edgeCompar1)),ncol(expectNet1Edges))
    expect_equal(sum(as.numeric(edgeCompar2)),ncol(expectNet2Edges))
    expect_equal(length(outputGraphs),2)
})

