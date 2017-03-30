#DIFFUSE Function

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

        #Check Class Mode and Type of AGG

        classModeTypeAGG <- c(class(AGG), mode(AGG), typeof(AGG))

        if (classModeTypeAGG != c("igraph", "list", "list")) {
            stop("input AGG is not an igraph object")
        }

        #Check other argmuments
        check1 <- .checkArgs(algorithm,character(length = 1),checkSize = TRUE)
        check2 <- .checkArgs(silent, logical(length = 1), checkSize = TRUE)
        check3 <- .checkArgs(noLog, logical(length = 1), checkSize = TRUE)
        check4 <- .checkArgs(param, list())

        if (c(check1,check2,check3,check4) != character(length = 0) ) {
            stop("Error: incorrect argument(s) for algorithm, silent, noLog, and/or param.
                 See documentation for DIFFUSE function")
        }

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
                        UUID = "12345",
                        input = "AGG",
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
