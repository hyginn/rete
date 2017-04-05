#DIFFUSE Function

#' Diffuse heat across edges and weight edges by heat
#'
#' @param AGG annotated gene graph (AGG), an igraph object with the following: vertices, with
#' gene name attribute, gene score attribute; edges (directed) with gene names corresponding to
#' 'source' and 'receiver' vertices
#' @param algorithm Method of distributing influence. Currently only "Leis" method, describing heat
#' diffusion process from Leisseron et al. (2015) available. Future versions of this function
#' may have additional methods. Once influence is distributed from a gene based upon its initial
#' 'gene score' (see score function for details), edges directed from the gene are assigned
#' the 'influence' calculated for them. See details section for more information
#' @param param list of parameters for algorithm, defaults to a list containing 'rete.beta',
#' a globally set option (see documentation for options()). See details section for more information
#' @param silent if TRUE, suppress printing of process to console
#' @param noLog if TRUE, suppress logging of run and object metadata to rete's log
#'
#' @return  EGG (igraph object with content equivalent to AGG, but with an influence attribute
#' for edges)
#'
#' @details The purpose of the function is to take an AGG input, and, depending on the algorithm,
#' assign directed influences from one vertex to another based on network topology and a given
#' vertex's gene score. For the Leis method, the following procedure is followed: 1) calculate
#' the matrix W(i,j) as 1/deg(j) if there is a directional edge from j to i, 2) calculate matrix F
#' as Beta * inverse(I-(1-Beta)*W), where W is the transition matrix as described before, I is
#' the identity matrix of dimensions (number of vertices x number of vertices), and Beta is an
#' 'insulation' parameter, 3) calculate matrix E as F*D, where D is a diagonal matrix where every
#' element (i,i) is the 'heat' score/gene score for vertex(i), 4) from each E(i,j) that is not
#' equal to 0, assign the value E(i,j) to the directional edge going from vertex j to vertex i.
#' A more precise definition of the model E=Beta*inverse(I-(1-Beta)*W)*D is a random walk where,
#' if one starts from node j, the probability of moving onto another given node will be
#' (1-Beta)*1/deg(j), where Beta is the probability of restarting from the start node. If this
#' process is repeated ad infinitum, a stationary distribution of F=Beta*inverse(I-(1-Beta)*W)
#' is reached, where in each column j of matrix F, a vector of probabilities is given for landing
#' on each element (i) in the column; i.e. element i of column j represents the probability of
#' reaching node i if a random walk starts from node j. Thus, element E(i,j) represents the
#' probability of reaching node i from node j when at equilibrium multiplied by the initial
#' 'heat' (gene score) put on node j. For all E(i,j) where there exists a directed edge from
#' j to i, edge j -> i is equal to the value of E(i,j).
#'
#' @export

DIFFUSE <- function(AGG = NULL, algorithm = "Leis", param = list(getOption("rete.beta")),
                    silent = FALSE, noLog = FALSE) {
        startTime <- Sys.time()


        if (!silent) {
            consoleVect <- paste("Diffuse started at", Sys.time(), sep = "")
            print(consoleVect) #print beginning of function to console
        }
#========= Checking Arguments =============================

        if (!silent) {
            consoleVect <- "Checking Arguments"
            print(consoleVect) #print beginning of function to console
        }

        #Check if AGG is an igraph object, and check if it's directed
        if (!(igraph::is.igraph(AGG) && igraph::is.directed(AGG))) {
            stop("Error: AGG must be an iGraph object with directed edges")
        }
        #Check if AGG is in fact an AGG object
        if (attr(AGG, "type") != "AGG") {
            stop("Error: input for AGG is not of type AGG")
        }

        #Check other argmuments
        check1 <- .checkArgs(algorithm,character(length = 1),checkSize = TRUE)
        check2 <- .checkArgs(silent, logical(length = 1), checkSize = TRUE)
        check3 <- .checkArgs(noLog, logical(length = 1), checkSize = TRUE)
        check4 <- .checkArgs(param, list())

        if (length(c(check1,check2,check3,check4)) > 0 ) {
            stop("Error: incorrect argument(s) for algorithm, silent, noLog, and/or param.
                 See documentation for DIFFUSE function")
        }


           #Assign AGG vertices and edges to data frames
           AGGverts<- igraph::as_data_frame(AGG, "vertices")
           AGGedges<- igraph::as_data_frame(AGG, "edges")

           sizeEdges <- nrow(AGGedges) #Number of edges, used for indexing purposes
           sizeVerts <- nrow(AGGverts) #Number of vertices, used for indexing purposes

    if (algorithm == "Leis") {


        if (!silent) {
            consoleVect<- "Method Selected: Leis"
            print(consoleVect)
        }

        #check beta
        if ( length(.checkArgs(param[[1]], numeric(length = 1))) > 0) {
            stop("Leis Method selected. Beta must be class numeric of type
                 double")
        }

        #creating transition matrix W
        HGNCsymb<-AGGverts$name #get a vector of AGG vertex names

        W<-matrix(data = numeric( length( sizeVerts )^2 ), nrow = sizeVerts , ncol = sizeVerts,
                  dimnames = list(HGNCsymb, HGNCsymb))

        vertDegree <- numeric(length = sizeVerts) #Will store degree of each vertex
        names(vertDegree) <- AGGverts$name


        if (!silent) {
            consoleVect<-"Calculating degree for vertices"
            print(consoleVect)
        }
        #calculating degree via following procedure



        for (i in 1:sizeVerts) {


            edgeIndex <- as.logical(grepl(HGNCsymb[i],AGGedges$from, fixed = TRUE)
                                    | grepl(HGNCsymb[i],AGGedges$to, fixed = TRUE))
            #generate a character vector of genes present in edges including vertex of interest
            connectedGenes <- c(AGGedges$from[edgeIndex],
                             AGGedges$to[edgeIndex])
            #reduce to unique elements
            connectedGenes <- unique(connectedGenes)
            #remove vertex of interest from list
            connectedGenes <- connectedGenes[-grep(HGNCsymb[i],connectedGenes, fixed = TRUE)]

            deg <- length(connectedGenes) #degree of vertex of interest

            vertDegree[grep(HGNCsymb[i],names(vertDegree), fixed = TRUE)] <- deg

        }

    #Assigning values to matrix W: See Leiseron et. al 2015 for procedure description


        if (!silent) {
            consoleVect <-  "Constructing matrix W"
            print(consoleVect)
        }


        for (n in 1:sizeVerts) {


           #look for directional edges with gene j as the source
           jEdgeIndex <- n

           #List of receivers of directed edge from j
           jTargetList <- AGGedges$to[jEdgeIndex]


               #W [i (receiver) , j (source)] is equal to 1/degree(j)
               W[jTargetList,n] <- 1/vertDegree[names(vertDegree) == jElement]




        } #end construction of W

    #======== Part 3 of Leis: calculating the matrix F ==============
        betaVal <- param[[1]] #set beta Value, or the 'insulation' coefficient

        #Following equation F= B*inv(I-(1-B)*W)

        consoleVect <- "calculating matrix F"
        if (!silent) {
            print(consoleVect)
        }

        idMat <- diag(nrow = sizeVerts, ncol = sizeVerts)
        preInvert <- idMat-(1-betaVal)*W #perform initial operations
        postInvert <- solve(preInvert) #invert matrix
        matrixF <- betaVal*postInvert #Finish and make matrix F

    #======== Part 4 of Leis: calculating matrix E ===================

        consoleVect <- "calculating matrix E"


        #Make a diagonal matrix D, or vector of heat scores in diagonal form

        D <- diag(x = AGGverts$Gene_Score)

        #Calculate matrix E (weighted diffusion matrix) as F*D (F postmultiplied
        #by D)

        matrixE <- matrixF %*% D

    #======== Part 5 of Leis: extracting 'heat' influence from matrix E =====

        consoleVect <- c(consoleVect, "Assigning heat/influence scores from E to edges")


        influenceVect <- numeric(length = sizeEdges)
        #vector of zeros, to be filled with E(i,j) values

        for (k in 1:sizeEdges) {
            j <- AGGedges$from[k]
            i <- AGGedges$to[k]

            influenceVect[k] <- matrixE[i,j]
        }

        AGGedges$Influence <- influenceVect

        EGG <- igraph::graph_from_data_frame(AGGedges, directed = TRUE,
                                             vertices = AGGverts)

    }
    else {
        stop("Error: No method selected")
    }

           finish <- Sys.time()

        #setup metadata
            attr(EGG , "type") <- "EGG"
            attr(EGG , "version") <- "1.0"
            attr(EGG , "UUID") <- uuid::UUIDgenerate()


#==================write the log=====================

    if (!noLog) {

        if (!silent) {
            consoleVect <- "Writing Log"
            print(consoleVect)
        }

        #write out the function call
        myCall <- character()
        myCall[1] <- "DIFFUSE("
        # ToDo ... update
        myCall[2] <- "AGG = AGG" #Don't know what to set this to
        myCall[3] <- paste("algorithm =", algorithm, sep = " ")
        myCall[3] <- paste("param =", param, sep = " ")
        myCall[4] <- paste("silent =", silent, sep = " ")
        myCall[5] <- paste("writeLog =", writeLog, sep = " ")

        myNotes <- character()

        myNotes <- c(myNotes, sprintf("AGG UUID: %s", attr(AGG, "UUID")))

        myNotes <- c(myNotes, sprintf("EGG UUID: %s", attr(EGG, "UUID")))

        myNotes <- c(myNotes, paste("start time", startTime, sep = " "))

        myNotes <- c(myNotes, paste("end time", finish, sep = " "))

        logEvent( eventTitle = "DIFFUSE",
                  eventCall = myCall,
                  input = character(), #put object name of AGG here?
                  notes = myNotes,
                  output = c("EGG"))

    }


          if (!silent) {
              consoleVect <- "Finished, returning EGG"
              print(consoleVect)

          }

           return(EGG)



}
#[END]
