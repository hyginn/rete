FINDSUBv1_2 <- function(MethType="Leis",Thresh,EGG,minOrd,logResults=FALSE,silent=FALSE) {
#NOTE: THIS IS THE MOST RECENT FORM OF FINDSUB, AND HAS BEEN TESTED AND SHOWN TO WORK IN SCRIPTS
#TITLED FINSUB_Test1.R through FINDSUB


    #INPUTS:
    #   MethType: method of edge removal
    #   Thresh: threshold for edge removal for given MethType
    #   EGG: 'Heat Diffusion' annotated AGG, iGraph object.
    #   minOrd: minimum number of vertices in subnetwork for it to be included
    #   logResults: logical. do we log our results to an Rdata file?
    #   silent: logical. do we print the progress of the function to console?
    #OUTPUT:
    #   graphList: List of subnetworks resulting after edge removal and
    #       subnetwork construction

#=====================================================================================
#Module 1: setting vertex and edge dataframes, performing checks

    logVect<-"Begin function FINDSUB" #begin log
    messt<-paste("Begun at date/time ",as.character(Sys.time()))
    logVect<-c(logVect,messt) #catenate log messages

    if (silent==FALSE) {
        print(logVect[length(logVect)])
    }
    logName<-paste("FINDSUB_Run_",as.character(Sys.Date()),".Rdata",sep="") #Set a name for the log

    #Now all the Checks

    #Have we downloaded iGraph
    if(!as.logical(require(igraph))) { #you need iGraph functions for this function.
        stop("Error: Need igraph package.")
    }
    #Is minOrd of correct type class and mode?
    if ((typeof(minOrd)=="integer" || typeof(minOrd)=="double") && (as.integer(minOrd)==minOrd)) {
        if (minOrd<2) {
            stop("minOrd must be 2 or greater")
        }
    }
    else {
        stop("minOrd must be an integer number, of type integer or double")
    }
    #is EGG an iGraph object?
    if (class(EGG) != "igraph") {
        stop("Error: EGG must be an iGraph object")
    }
    #Do we have a method type? For the time being we only have Leis
    if (MethType!="Leis") {
        stop("Error: MethType for edge removal
             must be Leis until later version add new method types")
    }
    #is thresh numeric and of class double, and a postive number (or at least 0)
    if (typeof(Thresh)== "double" && mode(Thresh) =="numeric") {
        if (Thresh<0) {
            stop("Thresh must be greater than 0!")
        }
    }
    else {
        stop("Thresh must be of type double and mode and class numeric")
    }
    #Is logResults a logical?
    if (class(logResults)!="logical") {
        stop("Error: logResults must be of mode type and class logical")
    }
    #is silent logical
    if (class(silent)!="logical") {
        stop("Error: silent must be of mode type and class logical")
    }

    vertEGG<-as_data_frame(EGG,what="vertices")
    names(vertEGG)[1][names(vertEGG)[1]=="name"]<-"HGNC_symbol"#renaming a column
    edgeEGG<-as_data_frame(EGG,what="edges")

    #inputs<-data.frame(MethType,Thresh,vertEGG,edgeEGG,minOrd,logResults,silent)
    #inputFname<-paste("Inputs,Run_",as.character(Sys.Date()),".RData")
    #save(inputs,inputFname)

    #updating log, below
    mess1<-paste("vertEGG assigned: No. Vertices=",as.character(nrow(vertEGG)))
    logVect<-c(logVect,mess1)
    if (silent==FALSE) {
        print(logVect[length(logVect)])
    }
    mess2<-paste("edgeEGG assigned: No. Edges= ",as.character(nrow(edgeEGG)))
    logVect<-c(logVect,mess2)
    if (silent==FALSE) {
        print(logVect[length(logVect)])
    }


#=====================================================================================
#Module 2: Method Selection and Edge Removal

    mess3<-paste("Method Selected (MethType:",MethType)
    logVect<-c(logVect,"========Begin Module 2: Edge Removal========",mess3)
    if (silent==FALSE) {
        print(logVect[length(logVect)-1])
        print(logVect[length(logVect)])
    }

    if (MethType == "Leis") {
    #=================SubModule: Leis Method of edge removal==========================

            removalList<-vector(mode="logical",length=nrow(edgeEGG))
            #Approach: remove
            for (i in 1:nrow(edgeEGG)) {
            #Need method of indexing edges
                messL1<-paste("Edge Removal, Iteration ",as.character(i))
                logVect<-c(logVect,messL1)
                if (silent==FALSE) {
                    print(logVect[length(logVect)])
                }

            edgeScore<-edgeEGG[i,"Influence"] #assign edge score of directed edge N
                if (edgeScore<Thresh) {
                    #if edgeScore less than threshold, put index into removalList
                    #as TRUE. removalList will then
                    removalList[i]=TRUE
                }

            }
            #indices of sub threshold edges found, now remove them
            removedEdges<-edgeEGG[removalList==TRUE,c("edgeID","from","to")]
            #removedEdges will be displayed
            noRemoved<-sum(as.numeric(removalList))
            edgeEGG<-edgeEGG[removalList==FALSE,]#actual step that removes edges
            messL2<-paste("Edge Removal Complete:",as.character(noRemoved),"edges removed")
            logVect<-c(logVect,messL2)
            if (silent==FALSE) {
                print(logVect[length(logVect)])
                print(removedEdges)
            }

            vertRemove<-vector(mode="logical",length=nrow(vertEGG))
            for (i in 1:nrow(vertEGG)) {
                protName<-as.character(vertEGG[i,"HGNC_symbol"]) #name refers to HGNC symbol
                if (sum(as.numeric(grepl(protName,edgeEGG$from)))
                    +sum(as.numeric(grepl(protName,edgeEGG$to))) == 0) {
                    vertRemove[i]=TRUE
                }

            }
            removedVertices<-vertEGG[vertRemove==TRUE,'HGNC_symbol']
            noRemovedVertices<-sum(as.numeric(vertRemove))
            messL3<-paste("Vertex Removeal Complete:", as.character(noRemovedVertices),"vertices removed")
            logVect<-c(logVect,messL3)
            if (silent==FALSE) {
                print(logVect[length(logVect)])
                print(removedVertices)
            }

            vertEGG=vertEGG[vertRemove==FALSE,]  #remove all vertices for which
                                                # no edges found

    }
    else {
        #update log and write it
        logVect=c(logVect,"Error: no method selected, ending function")
        logFINDSUB<-data.frame(LogFINDSUB=logVect,stringsAsFactors = FALSE)
        if (silent==FALSE) {
            print(logVect[length(logVect)])
        }
        save(logFINDSUB,file=logName)#This will have to be edited to save to a differnt location, perhaps
        stop("No Method selected for MethType")
    }




#=======================================================================================
#Module 3: subnetwork contruction

#idea: iteratively add edges to a subnetwork until no more edges with at least 1 mutal
    # vertex to one of teh subnetwork edges are found from the edgeEGG subset
    # repeat process with first of remaining edges

    #In constructing subnetworks, start out with an edge serving as a seed (must not have been previously
    #used as a seed, or included in any existing subnetworks)
    #Then, from this seed, search for edges in "query" pool (not yet incorporated in subnetwork)
    #sharing at least one vertex with the seed
    #These edges are added to the subnetwork, and for each of the newly added edges
    #the search for edges sharing vertices is repeated in a pool of edges
    #that have not been incorporated into the subnetwork yet
    #This process is repeated as long as new edges are discovered and aded to the subnetwork
    #When no new edges are 'discovered' this way, the subnetwork is complete

    mess4<-"==============Begin Module 3: Subnetwork Construction======================"
    logVect<-c(logVect,mess4)
    if (silent==FALSE) {
        print(logVect[length(logVect)])
    }
    #Approach of this module: go through the dataframe of edges edgesEGG. Take edge from first row
    #then look amongst a pool containing all other edges for edges that share at least one vertex.
    #If you find such edges, they become 'newly added' edges to your subnetwork,
    #and for each of these you repeat the
    #process of searching in a query pool of edges that have not been used as a 'seed edge (explanation
    #will follow)' and adding edges where a mutual vertex exists. You repeat this step however many
    #times it takes before you are unable to find additional edges. It should be noted that each time
    #a 'query' is added to the subnetwork, its index in edgesEGG is taken, and it is excluded from
    #the query pool in the next search among queries(i.e. the next time you try to match one of the
    #newly added edges to one of the query edges)

    graphList<-list() #This will store the igraph objects



    includedList<-vector(mode="logical",length=nrow(edgeEGG)) #Indices of edges that
                                                             #have made it into a subnetwork, regardless
                                                                #of whether or not the network meets
                                                             #minimum order threshold


    seedList<-vector(mode="logical",length=nrow(edgeEGG))
            #seedList= vector of indices, TRUE if indexed edge used as a seed. If used as a seed, either it's
            #included or will not find any partner edges. This vector exists so that we have a list of
            #things that were used as seeds but might not have been incorporated into a network of
            #appropriate size. If you haven't been incorporated into a network of appropriate size as a
            #seed, by definition, you've constructed the largest possible network containing that edge
            # and it still hasn't gotten past the threshold network size


    for (i in 1:nrow(edgeEGG)) {
        mess4.1<-paste("Iteration ",as.character(i))
        logVect<-c(logVect,mess4.1)
        if (silent==FALSE) {
            print(logVect[length(logVect)])
        }
        if (includedList[i]==FALSE) { #If indexed edge is not yet in a subnetwork



            seedEdgeID<-edgeEGG[i,"edgeID"] #retrieve identifier so we remove seed edge from query pool
            edgePartners<-edgeEGG[i,c("from","to")]
            mess4.2<-paste("Seeding Edge of edgeID ",seedEdgeID,"with interactors",
                           as.character(edgePartners[1])," & ",as.character(edgePartners[2]))
            logVect<-c(logVect,mess4.2)
            if (silent==FALSE) {
                print(logVect[length(logVect)])
            }


            subNetEdges<-edgeEGG[i,] #This sets up a variable containing
                                     #1st edge of new subnetwork (i.e. seed)
            seedList[i]<-TRUE

            cont<-TRUE #logical: do we continue below while loop
            #Basically, with the first edge, we search for anything with a mutual vertex
            #If we find anything that satisfies that condition, we add it to the subnetwork
            #cont=TRUE in this case
            #If we don't find anything that satisfies that condition, cont=FALSE
            #If we added a new edge, we repeat the loop, and see if any edges have mutal vertices
            #We repeat this process until no new edges are found


            subNetTried<-vector(length=1) # used to distinguish just added edges from edges added in
                                          #at a previous iteration

            M=0 #keeps track of how many queries have yielded an edge with an overlapping vertex

            while (cont==TRUE) {
                # Here, we generate a list of 'new edges' from which to search for other edges
                # that share mutual vertex. I.e. if we haven't tried an edge in the subnetwork
                # to find edges sharing a mutual vertex, we try to do just that, and once
                # we have tried that for that subnetwork edge, we add its index to the
                #vector subNetTried, and we don't use it again to search for connected edges
                N<-0
                newEdges<-subNetEdges[subNetTried==FALSE,] #only subnet edges that have not been tried
                edges2Add<-data.frame()

                for (n in 1:nrow(newEdges)) {
                    InEdgeVertices<-newEdges[n,c("from","to")] #changed 3.1.17. from subNetEdges to newEdges
                    queries<-edgeEGG[includedList==FALSE & seedList==FALSE,] #included list and seedlist
                    #are logical vectors that get updated to included edges that get included in subnetworks
                    # and edges that are used as seeds, respectively
                    #edited 2.28.17
                    if (silent==FALSE) {
                        print("Remaining Queries")
                        print(queries)
                    }

                    if (nrow(queries) != 0){
                        for (m in 1:nrow(queries)) { #edited 2.28.17

                            queryVertices<-queries[m,c("from","to")]
                            #If we see any pair of vertices sharing same gene name, we get a logical
                            #true, which converts to a numeric 1
                            cond1<-InEdgeVertices[1]==queryVertices[1]
                            cond2<-InEdgeVertices[1]==queryVertices[2]
                            cond3<-InEdgeVertices[2]==queryVertices[1]
                            cond4<-InEdgeVertices[2]==queryVertices[2]
                            vectCond=c(cond1,cond2,cond3,cond4)
                            sumVectCond<-sum(as.numeric(vectCond))
                            sumVectCondLogical<-(sumVectCond>0)
                            if (sumVectCondLogical==TRUE) {
                                #If any of the cond variables are
                                #true, then the sum of TRUES and FALsES expressed as 1s and 0s
                                #should be greater than 0
                                #In this case, add the edge to the variable "edges2Add"
                                #which contains the edges to add upon completion of the search
                                #for edges connected to the newEdges (edges added on previous iteration)
                                #or the 'seed' edge of the subnetwork if this is the first iteration of the
                                #while loop
                                edges2Add<-rbind(edges2Add,queries[m,])

                                #We also get the edge ID and add it to our index of included edges
                                #so we don't try to 'seed' a new network from something that
                                #already exists in one

                                includedEdgeID<-queries[m,"edgeID"]
                                includedListIndex<-edgeEGG[,"edgeID"]==queries[m,"edgeID"]
                                includedList[includedListIndex]=TRUE #indices by which to change
                                #query data frame indicated (edges to include from edgeEGG)
                                #but new indices not appied until the next iteraction of n

                                # We also want to remove the current query from unused queries


                                M=M+1
                                N<-N+1
                            } #end of if (sumVectCondLogical) (Do we have any mutual vertices?)
                        }#end of for m=1:nrow queries
                    } #end of if (nrows queries !=0 )

                } #end of for n=1:nrow newEdges
                #Now the newEdges aren't so new anymore: they've been used!
                #Set all current indices of subNetTried to TRUE
                lenSNT<-length(subNetTried)
                subNetTried[1:lenSNT]<-TRUE
                #expand subNetEdges to include the new edges, edges2Add
                subNetEdges<-rbind(subNetEdges,edges2Add)

                #Append a vector of 'FALSE's to subNetTried
                subNetTried<-c(subNetTried,vector(length=nrow(edges2Add)))

                #reminder: subNetTried serves as an index to distinguish 'just added'
                #edges from those from the previous iteration of while loop

                mess4.3<-paste("On search iteration ",as.character(M)," in subnetwork for
                               seed ",
                               as.character(i)," found ",as.character(nrow(edges2Add)),
                               " new edges connected to subnetwork")
                logVect<-c(logVect,mess4.3)

                if (silent==FALSE) {
                    print(logVect[length(logVect)])
                    print("New edges")
                    print(edges2Add)
                }

                if (N==0) { #If you haven't found anything new from your new edges
                            #end the while loop
                    cont=FALSE
                }

            } #end of while loop






            #In this section of code, we set up a data frame of vertices of the subnetwork
            subNetVertices<-data.frame() #create empty data frame to store vertex data
            subNetVertexList<-c(as.character(subNetEdges[,"from"]),as.character(subNetEdges[,"to"])) #vector of vertex IDs
            subNetVertexList<-unique(as.character(subNetVertexList))#remove duplicate elements

             if (length(subNetVertexList)>=minOrd) {#if the number of vertices>=minimum Order for subnetworks
                mess4.3<-paste("Order of subNetwork from seed ",as.character(i),">=minimum Order.",
                               " creating data frame with appropriate vertices")
                logVect=c(logVect,mess4.3)

                includedList[i]=TRUE #Marking index i of EGG as an included edge (i.e. seed edge)

                if (silent==FALSE) {
                    print(logVect[length(logVect)])
                    print("subnetwork vertex HGNC symbols=")
                    print(subNetVertexList)
                }
                for (n in 1:length(subNetVertexList)) {#Go through the list of vertices (gene names)
                    for (m in 1:nrow(vertEGG)) {
                        if (vertEGG[m,"HGNC_symbol"]==subNetVertexList[n]) {
                            #if we find a name match in vertices for name at index n, add the vertex data to the data frame
                            subNetVertices<-rbind(subNetVertices,vertEGG[m,])
                        }
                    }

                }

                #Then make the iGraph object, and attach it to graphList
                g<-graph_from_data_frame(subNetEdges,directed=TRUE,vertices=subNetVertices)
                graphList<-c(graphList,list(g))

                mess4.4<-paste("iGraph object for seed ",as.character(i), " created, with",
                               as.character(nrow(subNetEdges))," edges and ",
                               as.character(nrow(subNetVertices)), " vertices")
                logVect<-c(logVect,mess4.4)
                if (silent==FALSE) {
                    print(logVect[length(logVect)])
                }

            }
            else {
                mess4.5<-paste("Network from seed ",as.character(i)," did not meet minOrd threshold.
                               Subnetwork not converted to iGraph object")
                logVect<-c(logVect,mess4.5)
                if (silent==FALSE) {
                    print(logVect[length(logVect)])
                }
            }#end of if (numVertices>= minOrd) else ...



        } #end of embedded if statement that starts the subnetwork seed

    } #end of for loop going through edgeEGG
    numNetworks<-length(graphList)
    mess4.5<-paste("function complete, with ", as.character(numNetworks)," subnetworks produced and stored in iGraph objects")
    logVect<-c(logVect,mess4.5)
    if (silent==FALSE) {
        print(logVect[length(logVect)])
    }

    if (logResults==TRUE) {
        logFINDSUB<-data.frame(LogFINDSUB=logVect)
        save(logFINDSUB,file=logName)

    }


    return(graphList)

} #end of function