#findSub Tests

testthat::context("Test1 for Network Discovery")

testthat::test_that("2 Networks form of expected size, have correct heat
                    score assigned", {

                        #Note that Lines 11:47 are repeated in each context so that
                        #variables are reassigned, and do not carry over from previous contexts

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
                        #MethType: method of edge removal; Thresh: influence threshold; inputgraph: iGraph object containing
                        #network with assigned 'influence' to each edge; logResults: do we save a log of the process? ; silent
                        #if false, print what appears in the process log while the process is in action, as well as additional
                        #material (e.g. list edges removed, vertices removed)

                        #below 4 lines manually constructing expected networks
                        expectNet1Edges<-BigNetEdges[1:10,]
                        expectNet1Edges<-expectNet1Edges[order(expectNet1Edges$edgeID),]
                        expectNet1Verts<-BigNetVertices[1:4,]
                        expectNet2Edges<-BigNetEdges[c(15:16,19:20),]
                        expectNet2Edges<-expectNet2Edges[order(expectNet2Edges$edgeID),]
                        expectNet2Verts<-BigNetVertices[5:7,]
                        #below 6 lines loading iGraph object contents from output into data frames
                        netwk1<-outputGraphs[[1]]
                        netwk1Edges<-igraph::as_data_frame(netwk1,what="edges")
                        netwk1Verts<-igraph::as_data_frame(netwk1,what="vertices")
                        netwk2<-outputGraphs[[2]]
                        netwk2Edges<-igraph::as_data_frame(netwk2,what="edges")
                        netwk2Verts<-igraph::as_data_frame(netwk2,what="vertices")

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
                        testthat::expect_equal(nrow(expectNet1Verts),nrow(netwk1Verts))
                        testthat::expect_equal(nrow(expectNet2Verts),nrow(netwk2Verts))
                        testthat::expect_equal(sum(as.numeric(vertCompar1)),nrow(expectNet1Verts))
                        testthat::expect_equal(sum(as.numeric(vertCompar2)),nrow(expectNet2Verts))
                        testthat::expect_equal(sum(as.numeric(edgeCompar1)),ncol(expectNet1Edges))
                        testthat::expect_equal(sum(as.numeric(edgeCompar2)),ncol(expectNet2Edges))
                        testthat::expect_equal(length(outputGraphs),2)
                        testthat::expect_equal(igraph::graph_attr(netwk1,"aggHeatScore"),22)
                        testthat::expect_equal(igraph::graph_attr(netwk2,"aggHeatScore"),17)
                    })
testthat::context("findSub Test 2: testing effect of raising minOrd to 4")
testthat::test_that("1 subnetwork formed, with appropriate aggregate heat score, edges, and vertices", {
    #Here we only expect 1 subnetwork to form, with ABC1, DEF2, GHI4, and JKL3

    #again, the setup of variables
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

    #end assignment of common variables

    #begin: changed inputs for test 2

    igraph::graph_attr(EGG,"delta")<-4 #setting delta parameter (score threshold for edge inclusion)

    minOrd <- 4

    #We should only get 1 subnetwork this time

    #run findSub
    outputGraphs<-findSub(method="Leis",EGG,minOrd,noLog = FALSE,silent=FALSE)




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
    netwk1Edges<-igraph::as_data_frame(netwk1,what="edges")
    netwk1Verts<-igraph::as_data_frame(netwk1,what="vertices")


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

    #checking: expected vertices and edges, expected number of networks formed
    testthat::expect_equal(nrow(expectNet1Verts),nrow(netwk1Verts))
    testthat::expect_equal(sum(as.numeric(vertCompar1)),nrow(expectNet1Verts))
    testthat::expect_equal(sum(as.numeric(edgeCompar1)),ncol(expectNet1Edges))
    testthat::expect_equal(length(outputGraphs),1)

    #checking: expected aggregate heat score
    testthat::expect_equal(igraph::graph_attr(netwk1,"aggHeatScore"),22)
})

testthat::context("checking the effect of lowered delta and minOrd")

testthat::test_that("3 subnetworks formed with appropriate edges, vertices, heat scores", {

    #Begin  setup of variables common to all tests
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

    #end setup of common variables

    #begin setup of other inputs for test

    #delta set to 3.1, minOrd to 2

    igraph::graph_attr(EGG,"delta")<-3.1 #setting delta parameter (score threshold for edge inclusion)

    minOrd <- 2

    #Here, we expect three subnetworks to form. UVW8 and XYZ9 have directed edges going each way of Influence 3.5
    #Thus, they were eliminated in Test1 with a delta of 4

    #run findSub
    outputGraphs<-findSub(method="Leis",EGG,minOrd,noLog = TRUE,silent=FALSE)
    #MethType: method of edge removal; Thresh: influence threshold; inputgraph: iGraph object containing
    #network with assigned 'influence' to each edge; logResults: do we save a log of the process? ; silent
    #if false, print what appears in the process log while the process is in action, as well as additional
    #material (e.g. list edges removed, vertices removed)

    #Note

    #With a threshold influence of 3.5, LMN5 OPQ6 RST7 should form a network,
    #ABC1,DEF2,GHI4,JKL3 should form a network
    #When this is lowered to 3, LMN5,OPQ6,RST7,UVW8,and XYZ9 should form a network
    #ABC1 DEF2 JKL3 GHI4 from same subnetwork
    #Lower further to 2, and the 'discovered subnetwork' will include all vertices, but not all edges
    #Lower influence threshold (delta) to 1, and every edge included

    #below 6 lines: set up expected edges and vertices
    expectNet1Edges<-BigNetEdges[1:10,]
    expectNet1Edges<-expectNet1Edges[order(expectNet1Edges$edgeID),]
    expectNet1Verts<-BigNetVertices[1:4,]
    expectNet2Edges<-BigNetEdges[c(15:16,19:20),]
    expectNet2Edges<-expectNet2Edges[order(expectNet2Edges$edgeID),]
    expectNet2Verts<-BigNetVertices[5:7,]
    expectNet3Edges<-BigNetEdges[21:22,]
    expectNet3Edges<-expectNet3Edges[order(expectNet3Edges$edgeID),]
    expectNet3Verts<-BigNetVertices[8:9,]

    #below 9 lines: set up edges and vertices from output as data frames
    netwk1<-outputGraphs[[1]]
    netwk1Edges<-igraph::as_data_frame(netwk1,what="edges")
    netwk1Verts<-igraph::as_data_frame(netwk1,what="vertices")
    netwk2<-outputGraphs[[2]]
    netwk2Edges<-igraph::as_data_frame(netwk2,what="edges")
    netwk2Verts<-igraph::as_data_frame(netwk2,what="vertices")
    netwk3<-outputGraphs[[3]]
    netwk3Edges<-igraph::as_data_frame(netwk3,what="edges")
    netwk3Verts<-igraph::as_data_frame(netwk3,what="vertices")

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

    edgeCompar3<-vector(length=ncol(expectNet3Edges))
    if (ncol(expectNet3Edges)==ncol(netwk3Edges) && nrow(expectNet3Edges)==nrow(netwk3Edges)) {
        for (i in 1:ncol(expectNet2Edges)) {
            truRows<-expectNet3Edges[,i]==netwk3Edges[,i]
            sumTruRows<-sum(as.numeric(truRows))
            if (sumTruRows==nrow(expectNet3Edges)) {
                edgeCompar3[i]<-TRUE
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
    vertCompar3<-vector(length=nrow(expectNet3Verts))
    for (i in 1:nrow(expectNet3Verts)) {
        ind3<-grepl(expectNet3Verts$HGNC_symbol[i],netwk3Verts$name)
        if (sum(as.numeric(ind1))>0) {
            if (netwk3Verts$Gene_Score[i]==expectNet3Verts$Gene_Score[ind3]) {
                vertCompar3[i]<-TRUE
            }
        }
    }

    #do we have 3 networks present with expected edges and vertices

    testthat::expect_equal(nrow(expectNet1Verts),nrow(netwk1Verts))
    testthat::expect_equal(nrow(expectNet2Verts),nrow(netwk2Verts))
    testthat::expect_equal(nrow(expectNet3Verts),nrow(netwk3Verts))
    testthat::expect_equal(sum(as.numeric(vertCompar1)),nrow(expectNet1Verts))
    testthat::expect_equal(sum(as.numeric(vertCompar2)),nrow(expectNet2Verts))
    testthat::expect_equal(sum(as.numeric(vertCompar3)),nrow(expectNet3Verts))
    testthat::expect_equal(sum(as.numeric(edgeCompar1)),ncol(expectNet1Edges))
    testthat::expect_equal(sum(as.numeric(edgeCompar2)),ncol(expectNet2Edges))
    testthat::expect_equal(sum(as.numeric(edgeCompar3)),ncol(expectNet3Edges))
    testthat::expect_equal(length(outputGraphs),3)

    #do we have the expected heat scores for each network
    testthat::expect_equal(igraph::graph_attr(netwk1,"aggHeatScore"),22)
    testthat::expect_equal(igraph::graph_attr(netwk2,"aggHeatScore"),17)
    testthat::expect_equal(igraph::graph_attr(netwk3,"aggHeatScore"),14)
})

testthat::context("Test 4, part 1: Do we fail to get a network created if not strongly connected")
testthat::test_that("1 subnetwork forms with expected edges and vertices and heat score", {
    #Begin construction of common variables
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

    #end construction of common variables

    #Here, we will be assigning an edge influence attribute of 5 to edge GHI4->XYZ9. Since this does not
    #create a 'strongly connected' network that liks XYZ9 & UVW8 to ABC1:GHI4,
    #ABC1 DEF2 JKL3 GHI4 should form a network of 4 nodes and LMN5, OPQ6, RST7 should form
    #a subnetwork of 3 nodes.


    #Influence is reassigned to several edges such that after removing edges below the delta threshold,
    #ABC1 through GHI4 and UVW8:XYZ9 are liked by 1 directional edge. This does not create a strongly
    #connected graph that includes all of those vertices
    Influence[(From=="UVW8")&(To=="XYZ9")]=5
    Influence[(From=="XYZ9")&(To=="UVW8")]=5
    Influence[(From=="GHI4")&(To=="XYZ9")]=5

    BigNetEdges<-data.frame(from=From,to=To,Influence=Influence,edgeID=edgeID) #remaking the data frame
    inputGraph<-igraph::graph_from_data_frame(BigNetEdges,directed=TRUE,vertices=BigNetVertices) #and inputGraph
    igraph::graph_attr(EGG,"delta")<-4 #setting delta parameter (score threshold for edge inclusion)

    minOrd <- 3  #setting minOrd to 3s

    #below, applying findSub function
    outputGraphs<-findSub(method="Leis",EGG,minOrd,noLog=FALSE,silent=FALSE)

    #MethType: method of edge removal; Thresh: influence threshold; inputgraph: iGraph object containing
    #network with assigned 'influence' to each edge; logResults: do we save a log of the process? ; silent
    #if false, print what appears in the process log while the process is in action, as well as additional
    #material (e.g. list edges removed, vertices removed)



    #below 6 lines: set up expected edges and vertices
    expectNet1Edges<-BigNetEdges[c(1:10),]
    expectNet1Edges<-expectNet1Edges[order(expectNet1Edges$edgeID),]
    expectNet1Verts<-BigNetVertices[c(1:4),]
    expectNet2Edges<-BigNetEdges[c(15:16,19:20),]
    exoectNet2Edges<-expectNet2Edges[order(expectNet2Edges$edgeID)]
    expectNet2Verts<-BigNetVertices[5:7,]
    #below 9 lines: set up edges and vertices from output as data frames
    netwk1<-outputGraphs[[1]]
    netwk1Edges<-igraph::as_data_frame(netwk1,what="edges")
    netwk1Verts<-igraph::as_data_frame(netwk1,what="vertices")
    netwk2<-outputGraphs[[2]]
    netwk2Edges<-igraph::as_data_frame(netwk2,what="edges")
    netwk2Verts<-igraph::as_data_frame(netwk2,what="vertices")


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
    #Have we got the expected subnetworks?
    testthat::expect_equal(nrow(expectNet1Verts),nrow(netwk1Verts))
    testthat::expect_equal(nrow(expectNet2Verts),nrow(netwk2Verts))
    testthat::expect_equal(sum(as.numeric(vertCompar1)),nrow(expectNet1Verts))
    testthat::expect_equal(sum(as.numeric(vertCompar2)),nrow(expectNet2Verts))
    testthat::expect_equal(sum(as.numeric(edgeCompar1)),ncol(expectNet1Edges))
    testthat::expect_equal(sum(as.numeric(edgeCompar2)),ncol(expectNet2Edges))
    testthat::expect_equal(length(outputGraphs),2)
    #and do they have the right aggregate heat scores?
    testthat::expect_equal(igraph::graph_attr(netwk1,"aggHeatScore"),22)
    testthat::expect_equal(igraph::graph_attr(netwk2,"aggHeatScore"),17)
})

testthat::context ("Test 4, Part 2: When we make the previously 'not strongly' connected network strongly connected,
         does it get constructed?")
testthat::test_that("2 large networks found and are appropriately formed and scored", {

             #once again setting up common variables
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

             #end setup of common variables

             #Like FINDSUB_Test4, but setting edge weights such that after edge removal, a path will exist between
             #GHI4 and XYZ9 going in both directions, and consequently UVW8 and XYZ9 are included in a strongly connected
             #subnetwork with ABC1:GHI4

             #Here, we will be assigning an edge influence attribute of 5 to edge GHI4->XYZ9, and an
             #edge influence attribute of 5 to edge XYZ0->GHI4
             #ABC1 DEF2 JKL3 GHI4 and  should form a network of 6 nodes and LMN5, OPQ6, RST7 should form
             #a subnetwork of 3 nodes.

             #Another of this test is to show that construction of subnetworks can involve grouping
             #edges that are not grouped together i.e. consecutive in the list of edges.

             Influence[(From=="UVW8")&(To=="XYZ9")]=5
             Influence[(From=="XYZ9")&(To=="UVW8")]=5
             Influence[(From=="GHI4")&(To=="XYZ9")]=5
             Influence[(From=="XYZ9")&(To=="GHI4")]=5

             BigNetEdges<-data.frame(from=From,to=To,Influence=Influence,edgeID=edgeID) #remaking the data frame
             EGG<-igraph::graph_from_data_frame(BigNetEdges,directed=TRUE,vertices=BigNetVertices) #and inputGraph
             igraph::graph_attr(EGG,"EGGversion")<-"1.0"

             igraph::graph_attr(EGG,"delta")<-4 #setting delta parameter (score threshold for edge inclusion)

             minOrd <- 3  #setting minOrd to 3s

             #below, applying findSub function
             outputGraphs<-findSub(method="Leis",EGG,minOrd,noLog=FALSE,silent=FALSE)

             #MethType: method of edge removal; Thresh: influence threshold; inputgraph: iGraph object containing
             #network with assigned 'influence' to each edge; logResults: do we save a log of the process? ; silent
             #if false, print what appears in the process log while the process is in action, as well as additional
             #material (e.g. list edges removed, vertices removed)



             #below 6 lines: set up expected edges and vertices
             expectNet1Edges<-BigNetEdges[c(1:10,12,14,21:22),]
             expectNet1Edges<-expectNet1Edges[order(expectNet1Edges$edgeID),]
             expectNet1Verts<-BigNetVertices[c(1:4,8:9),]
             expectNet2Edges<-BigNetEdges[c(15:16,19:20),]
             exoectNet2Edges<-expectNet2Edges[order(expectNet2Edges$edgeID)]
             expectNet2Verts<-BigNetVertices[5:7,]
             #below 9 lines: set up edges and vertices from output as data frames
             netwk1<-outputGraphs[[1]]
             netwk1Edges<-igraph::as_data_frame(netwk1,what="edges")
             netwk1Verts<-igraph::as_data_frame(netwk1,what="vertices")
             netwk2<-outputGraphs[[2]]
             netwk2Edges<-igraph::as_data_frame(netwk2,what="edges")
             netwk2Verts<-igraph::as_data_frame(netwk2,what="vertices")


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

             #Do we have our two expected networks, one of them of size 7, with appropriate
             #edges and vertices?

             testthat::expect_equal(nrow(expectNet1Verts),nrow(netwk1Verts))
             testthat::expect_equal(nrow(expectNet2Verts),nrow(netwk2Verts))
             testthat::expect_equal(sum(as.numeric(vertCompar1)),nrow(expectNet1Verts))
             testthat::expect_equal(sum(as.numeric(vertCompar2)),nrow(expectNet2Verts))
             testthat::expect_equal(sum(as.numeric(edgeCompar1)),ncol(expectNet1Edges))
             testthat::expect_equal(sum(as.numeric(edgeCompar2)),ncol(expectNet2Edges))
             testthat::expect_equal(length(outputGraphs),2)

             #Do the networks constructed have appropriate aggreagate heat scores assigned?

             testthat::expect_equal(igraph::graph_attr(netwk1,"aggHeatScore"),36)
             testthat::expect_equal(igraph::graph_attr(netwk2,"aggHeatScore"),17)


         })

testthat::context("Test 5: Checking all the checks and that incorrect input results in stoppage of function")

#setup EGG and other variables for test data one more time

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

#end common variable setup

#begin Test 5 tests

testthat::test_that("Checking the checks", {
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=NULL)) #silent is null
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=2)) #silent is numeric
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=1,
                                       noLog = TRUE,silent="TRUE")) #silent with character

    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = NULL,silent=FALSE)) #noLog is null
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = "NULL",silent=FALSE)) #noLog is Character
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = 43,silent=FALSE)) #noLog is numeric

    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=-3,
                                       noLog = TRUE,silent=TRUE))
    testthat::expect_error(outputGraphs<-findSub(method="Blargh",EGG,minOrd=3, #character method
                                       noLog = TRUE,silent=TRUE))

    testthat::expect_error(outputGraphs<-findSub(method=NULL,EGG,minOrd=3, #null method
                                       noLog = TRUE,silent=TRUE))

    testthat::expect_error(outputGraphs<-findSub(method=2,EGG,minOrd=3, #numeric method
                                       noLog = TRUE,silent=TRUE))

    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=NULL, #null minord
                                       noLog = TRUE,silent=TRUE))



    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=2.2,
                                       noLog = TRUE,silent=TRUE)) #minOrd with non-integer

    testthat::expect_error(outputGraphs<-findSub(method=2,EGG,minOrd="Pizza",
                                       noLog = TRUE,silent=TRUE)) #MinOrd w character
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=1,
                                       noLog = TRUE,silent=TRUE)) #minOrd<2

})

#An added test: checking the version number

testthat::context("Checking that you can check version number")

igraph::graph_attr(EGG,"EGGversion")<-"1.2"

testthat::test_that("Wrong version number -> Error message", {
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::graph_attr(EGG,"EGGversion")<-1

testthat::test_that("EGGversion of type double -> Error message", {
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::graph_attr(EGG,"EGGversion")<-as.integer(1)

testthat::test_that("EGGversion of type integer -> Error message", {
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::remove.graph.attribute(EGG,"EGGversion")

testthat::test_that("EGGversion nonexistent -> Error message", {
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

testthat::context("Have a valid delta value")

igraph::graph_attr(EGG,"delta")<-"puppies"

testthat::test_that("Delta a character -> Error message", {
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::graph_attr(EGG,"delta")<-as.integer(1)

testthat::test_that("Delta an integer -> Error message", {
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::graph_attr(EGG,"delta")<-as.integer(1)

testthat::test_that("Delta an integer -> Error message", {
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::graph_attr(EGG,"delta")<-as.integer(1)

testthat::test_that("Delta an integer -> Error message", {
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})



testthat::test_that("No Delta-> Error message", {
    igraph::remove.graph.attribute(EGG,"delta")
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

#finally, no EGG! (Expect an error)
testthat::test_that("No EGG -> error", {
    igraph::graph_attr(EGG,"delta")<-4
    testthat::expect_error(outputGraphs<-findSub(method="Leis",NULL,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
    testthat::expect_error(outputGraphs<-findSub(method="Leis","EGG",minOrd=3,
                                       noLog = TRUE,silent=TRUE))
    testthat::expect_error(outputGraphs<-findSub(method="Leis",0,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
    rm(EGG)
    testthat::expect_error(outputGraphs<-findSub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})
