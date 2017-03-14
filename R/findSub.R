findSub <- function(method="Leis",EGG,minOrd,noLog = FALSE,silent = FALSE) {
#NOTE: THIS IS THE MOST RECENT FORM OF FINDSUB, AND HAS BEEN TESTED AND SHOWN TO WORK IN SCRIPTS
#TITLED FINSUB_Test1.R through FINDSUB

    #'Find Highly Connected Subnetworks after Eliminating 'Weak' Edges,
    #'and calculating aggregate 'heat' scores for these subnetworks
    #'
    #'@param method method of edge removal. Currently only "Leis" (Leiseron et al. 2015) method of edge
    #'removal available. See details for more information
    #'@param EGG an igraph object of format AGG with score annotations
    #'@param minOrd minimum number of vertices per generated subgraph
    #'@param noLog if FALSE, save process to an .Rdata log file
    #'@param silent if FALSE, print process to console in real time
    #'@return A list of igraph objects, with "Influence" edge attributes,
    #'"Gene Score" vertex attributes, "EGGversion" graph Attribute,
    #'"Aggregate Heat score graph attribute
    #'@details Depending on method used, edges are removed from EGG based upon the assigned 'delta' attribute
    #'and their 'Influence' scores(see next paragraph). In the "Leis" method (Leiseron et al. 2015), edges with
    #'Influence scores below the value of 'delta' are eliminated. Then, highly connected components of
    #'the EGG graph made of the remaining edges are found using igraph::components(),
    #'and subgraphs generated for these components. The vertices in these subgraphs come with the Gene Scores
    #'(score can differ depending on method), which for each subgraph's vertices are summed to find an
    #'aggregate heat score for each individual graph. E.g., if a subgraph's vertices had scores of 5,8,11
    #' the aggregate heat score would be 24. The aggregate heat score is then assigned as an attribute
    #' aggHeatScore. Subgraphs are stored in a list of subgraphs sN, which is returned. sN is ordered from
    #' highest aggHeatScore to lowest aggHeatScore of the subgraph.
    #' All subgraphs in sN contain same attributes as EGG and have aggHeatScore (aggregate heat score)
    #' added as an attribute.
    #'
    #'EGG must have the following: vertices: genes with HGNC symbol
    #'and Gene Score ('heat', see documentation for score function for more details) as the minimal vertex
    #'attributes; edges: directed edges from 1 gene to the next, derived from protein protein interaction data.
    #'Edges' 'from' and 'to' columns (edge vertices) named with HGNC symbols. Minimum attributes for edges
    #'include 'Influence', a score weighting edges depending on weighting method (for our purposes this
    #'will likely ). EGG must also come with a delta attribute (integer or of type double)
    #'(assigned threshold for edge removal regardless of method used in generating edge scores) and an
    #'EGGversion, a graph attribute with a character vector of length 1. Unless further changes are made,
    #'EGGversion must be "1.0".
    #'




    #INPUTS:
    #   method: method of edge removal
    #   delta: threshold for edge removal for given method
    #   EGG: 'Heat Diffusion' annotated AGG, iGraph object.
    #   minOrd: minimum number of vertices in subnetwork for it to be included
    #   noLog: logical. do we log our results to an Rdata file?
    #   silent: logical. do we print the progress of the function to console?
    #OUTPUT:
    #   sN: List of subnetworks resulting after edge removal and
    #       subnetwork construction

#=====================================================================================
#Module 1: Checks

    logVect<-"Begin function FINDSUB" #begin log
    messt<-paste("Begun at date/time ",as.character(Sys.time()),collapse="\n")
    logVect<-c(logVect,messt) #catenate log messages

    if (!silent) {
        print(logVect[length(logVect)])
    }

    #setting log name
    dateTimeName<-strsplit(as.character(Sys.time()), split="[- ]")
    yr<-dateTimeName[[1]][1]
    mm<-dateTimeName[[1]][2]
    dd<-dateTimeName[[1]][3]
    tm<-strsplit(dateTimeName[[1]][4],split="[:]")
    hr<-tm[[1]][1]
    min<-tm[[1]][2]
    sec<-tm[[1]][3]
    logName<-paste("FINDSUBrun",yr,"yr",mm,"month",dd,"day",hr,"hr",min,"min",
                   sec,"sec",".Rdata",sep = "") #Set a name for the log

    #Now all the Checks


    #Is minOrd of correct type class and mode?
    if ((typeof(minOrd) == "integer" || typeof(minOrd) == "double") && (as.integer(minOrd) == minOrd)) {
        if (minOrd < 2) {
            stop("minOrd must be 2 or greater")
        }
    }
    else {
        stop("minOrd must be an integer number, of type integer or double")
    }
#removed: check for iGraph object type
    #Do we have a method type? For the time being we only have Leis
    if (method != "Leis") {
        stop("Error: method for edge removal
             must be Leis until later version adds new method types")
    }
    #is numeric and of class double, and a postive number (or at least 0)
    #Is noLog a logical?
    if (class(noLog) != "logical") {
        stop("Error: noLog must be of mode type and class logical")
    }
    #is silent logical
    if (class(silent) != "logical") {
        stop("Error: silent must be of mode type and class logical")
    }
    #is EGG of correct EGGversion?
    if (igraph::graph_attr(EGG, "EGGversion")  !=  "1.0") {
        stop("Error: EGGversion must be 1.0")
    }


#=====================================================================================
#Module 2: Method Selection and Edge Removal

    mess3<-paste("Method Selected (method:",method,collapse="\n")
    logVect<-c(logVect,"========Begin Module 2: Edge Removal========",mess3)
    if (!silent) {
        print(logVect[length(logVect)-1])
        print(logVect[length(logVect)])
    }

    if (method == "Leis") {
    #=================SubModule: Leis Method of edge removal==========================
        #not quite sure if we should put delta assignment and checks upstream
        delta<- igraph::graph_attr(EGG,"delta") #assigning delta
        #assuming that EGG's delta is assigned as an attribute

        #check that delta is valid input
        deltaCond1<-(class(delta) == "numeric" && typeof(delta) == "double" && mode(delta) == "numeric")
        deltaCond2<-(delta >= 0)

        if (!deltaCond1 || !deltaCond2) {
            stop("delta must be of mode and type numeric, and class double, and must be greater than
                 or equal to 0")
        }

        #delete all edges with influence < delta
        EGG <- igraph::delete_edges(EGG, igraph::E(EGG)[igraph::E(EGG)$Influence < delta])

    }
    else {
        #update log and write it
        logVect = c(logVect,"Error: no method selected, ending function")
        logFINDSUB<-data.frame(LogFINDSUB = logVect,stringsAsFactors = FALSE)
        if (!silent) {
            print(logVect[length(logVect)])
        }
        save(logFINDSUB,file = logName)#This will have to be edited to save to a differnt location, perhaps
        stop("No Method selected for edge removal")
    }




#=======================================================================================
#Module 3: subnetwork contruction



    mess4<-"==============Begin Module 3: Subnetwork Construction======================"
    logVect<-c(logVect,mess4)
    if (!silent) {
        print(logVect[length(logVect)])
    }
#Taking Dr. Steipe's suggestions
    sN<-list()
    cEGG <- igraph::components(EGG, mode = "strong") #remaining components of EGG
    aggHeatScoreVect<-numeric()
    for (i in 1:cEGG$no) { #for every component (identified by an integer btw 1 and 'no')
                            #perform below procedure to create iGraph objects for
                            #strongly connected graph components with order > minOrd
       if (cEGG$csize[i] >= minOrd) {
           vertList<-list(names(cEGG$membership[cEGG$membership == i]))#list of vertex names
                                                                    #where vertices are members
                                                                    # of cluster i

           newGraph<-igraph::induced_subgraph( EGG , vertList[[1]] )
           # creating new iGraph object from vertices in a strongly connected component
           # then ordering newGraph by edgeID
           newGraphEdges<-igraph::as_data_frame(newGraph,what = "edges")
           newGraphVerts<-igraph::as_data_frame(newGraph,what = "vertices")
           edgeIDVect<-newGraphEdges$edgeID
           edgeIDVect<-order(edgeIDVect)
           newGraphEdges<-newGraphEdges[edgeIDVect,]
           newGraph<-igraph::graph_from_data_frame(newGraphEdges,vertices = newGraphVerts)


            #then calculating aggregate Heat Score of the included vertices
           aggHeatScore<-sum(igraph::vertex_attr(newGraph,"Gene_Score"))
           #and then assigning the aggHeatScore as an attribute of the graph
           igraph::graph_attr(newGraph,"aggHeatScore")<-aggHeatScore
           #and adding the value to a vector so we can reorder the sN upon
           #creation of all strongly connected subgraphs
           aggHeatScoreVect<-c(aggHeatScoreVect,aggHeatScore)

           sN<-c(sN,list(newGraph))#adding the graph to the sN
           mess4<-paste("Component", i ,"out of", cEGG$no , "is of order >= minOrd",collapse="\n")
           mess5<-paste("Created subgraph for component" , i , "and calculated aggregate heat score",
                        collapse="\n")
           logVect<-c(logVect,mess4,mess5)
           if (!silent) {
               print(logVect[length(logVect) - 1])
               print(logVect[length(logVect)])
           }


       }
    }
    #order subgraphs by aggregate heat score
graphOrder <- order(aggHeatScoreVect, decreasing = TRUE )
sN <- sN[graphOrder]



mess6<-paste("Function finished.", length(sN) , "strongly connected subgraphs found",
             collapse="\n")

logVect<-c(logVect,mess6)
if (!silent) {
    print(logVect[length(logVect)])
}

#if (!noLog) { #save log data if we choose noLog = FALSE
  #  save(logVect, file = logName)
#}
return(sN)

} #end of function

#[END]
