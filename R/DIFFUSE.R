#DIFFUSE Function

#' Diffuse heat across edges and weight edges by heat
#'
#' @param AGG annotated gene graph (AGG), an igraph object with the following: vertices, with
#' gene name attribute, gene score attribute; edges (directed) with gene names corresponding to
#' 'source' and 'reciever' vertices
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
#' process is repeated ad infinitum, an equilibrium distribution of F=Beta*inverse(I-(1-Beta)*W)
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

        consoleVect <- paste("Diffuse started at", Sys.time(), sep = "")
        if (!silent) {
            print(consoleVect) #print beginning of function to console
        }
#========= Checking Arguments =============================
        consoleVect <- "Checking Arguments"
        if (!silent) {
            print(consoleVect) #print beginning of function to console
        }

        #Check if AGG is an igraph object, and check if it's directed
        if (!(igraph::is.igraph(AGG) && igraph::is.directed(AGG))) {
            stop("Error: AGG must be an iGraph object with directed edges")
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

        #Extract AGG metadata for writing to log
        metaAGG <- attr(AGG,"meta")
           #Assign AGG vertices and edges to data frames
           AGGverts<- igraph::as_data_frame(AGG, "vertices")
           AGGedges<- igraph::as_data_frame(AGG, "edges")

           sizeEdges <- nrow(AGGedges) #Number of edges, used for indexing purposes
           sizeVerts <- nrow(AGGverts) #Number of vertices, used for indexing purposes

    if (algorithm == "Leis") {

        consoleVect<- "Method Selected: Leis"
        if (!silent) {
            print(consoleVect)
        }

        #check beta
        if ( length(.checkArgs(param[[1]], numeric(length = 1))) > 1) {
            stop("Leis Method selected. Beta must be class numeric of type
                 double")
        }

        #creating transition matrix W
        zeroVect <- numeric(length(sizeVerts)^2) #Make a vector of zeros from which to make W
        W<-matrix(data = zeroVect, nrow = sizeVerts , ncol = sizeVerts,
                  dimnames = list(AGGverts$name, AGGverts$name))

        vertDegree <- numeric(length = sizeVerts) #Will store degree of each vertex
        names(vertDegree) <- AGGverts$name

        consoleVect<-"Calculating degree for vertices"
        if (!silent) {
            print(consoleVect)
        }
        #calculating degree via following procedure

        HGNCsymb<-AGGverts$name #get a vector of AGG vertex names

        for (i in 1:sizeVerts) {

            consoleVect<-paste("vertex",i, "of", sizeVerts, sep = " ")
            if (!silent) {
                print(consoleVect)
            }


            edgeIndex <- as.logical(grepl(HGNCsymb[i],AGGedges$from, fixed = TRUE)
                                    | grepl(HGNCsymb[i],AGGedges$to, fixed = TRUE))
            #generate a list of genes present in edges including vertex of interest
            geneList <- c(AGGedges$from[edgeIndex],
                             AGGedges$to[edgeIndex])
            #reduce to unique elements
            geneList <- unique(geneList)
            #remove vertex of interest from list
            geneList <- geneList[-grep(HGNCsymb[i],geneList, fixed = TRUE)]

            deg <- length(geneList) #degree of vertex of interest

            vertDegree[grep(HGNCsymb[i],names(vertDegree), fixed = TRUE)] <- deg

        }

    #Assigning values to matrix W: See Leiseron et. al 2015 for procedure description

        consoleVect <-  "Constructing matrix W"
        if (!silent) {
            print(consoleVect)
        }

        doneEdges <- 0

        for (n in 1:sizeVerts) {



           jElement <- AGGverts$name[n]

           #look for directional edges with gene j as the source
           jEdgeIndex <- as.logical(grepl(jElement, AGGedges$from, fixed = TRUE))

           #List of receivers of directed edge from j
           jTargetList <- AGGedges$to[jEdgeIndex]

           for (m in 1:length(jTargetList)) {
               iElement <- jTargetList[m] #retrieve name of target of node j

               #W [i (receiver) , j (source)] is equal to 1/degree(j)
               W[iElement,jElement] <- 1/vertDegree[names(vertDegree) == jElement]

               doneEdges <- doneEdges + 1

               consoleVect <- paste("Have value 1/deg(j) for",
                                                   doneEdges,
                                                   "of", sizeEdges,
                                                   "assigned to W", sep = " ")
               if (!silent) {
                   print(consoleVect)
               }

           }

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

        consoleVect <- c(consoleVect, "calculating matrix E")


        #Make a diagonal matrix D, or vector of heat scores in diagonal form

        D <- diag(x = AGGverts$Gene_Score)
        dimnames(D) <- list(AGGverts$name,AGGverts$name)

        #Calculate matrix E (weighted diffusion matrix) as F*D (F postmultiplied
        #by D)

        matrixE <- matrixF%*%D

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
           meta <- list(version = "EGG_Version_1.0",
                        UUID = uuid::UUIDgenerate(),
                        input = paste("AGGuuID",metaAGG["UUID"][[1]], sep = ""),
                        time = finish)
          attr(EGG,"meta") <- meta
#==================write the log=====================

    if (!noLog) {
        consoleVect <- "Writing Log"
        if (!silent) {
            print(consoleVect[length(consoleVect)])
        }

        logVect <- paste("DIFFUSE started at", startTime, sep = "")
        logVect <- c(logVect, "Finished at", finish, sep = "")
        logVect <- c(logVect,paste("returned EGG object with", sizeEdges, "edges and",
                                   sizeVerts, "vertices", sep = " "))
        logVect <- c(logVect,"UUID and attributes",as.character(meta))
        logMessage(logVect)
    }

          consoleVect <- "Finished, returning EGG"
          if (!silent) {
              print(consoleVect[length(consoleVect)])
          }

           return(EGG)



}
#[END]
