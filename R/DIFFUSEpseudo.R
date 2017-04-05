# function (AGG, algorithm = "Leis", param = list(getOptions(rete.beta),)) , silent= FALSE
# noLog = FALSE) {
#
#
#
#       *Do all checks to make sure arguments are of right type
#
#
#       *If all checks pass, and if method = Leis {
#
#
#       AGGverts<- igraph::as_data_frame(AGG, "vertices")
#       AGGedges<- igraph::as_data_frame(AGG, "edges")
#
#       ===================#Now set up the matrix W=====================
#
#       sizeEdges <- nrow(AGGedges)
#       sizeVerts <- nrow(AGGverts)
#       zeroVect <- numeric(length(sizeVerts)^2) #Make a vector of zeros from which to make W
#       W<-matrix(data = zeroVect, nrow = sizeVerts , ncol = sizeVerts)
#
#       #assiging values to elements of W
#
#           #Prior to for loop, get degree of all vertices
#           degreeVect<-igraph::degree(AGG vertices)   #return numeric vector of vertex degrees
#
#       for ( n in 1:sizeEdges) {
#
#
#           iElement=AGGedges[n,from]
#           jElement=AGGedges[n,to]
#
#           #Find iElement and jElement in AGGverts to obtain their indices
#           iIndex <- AGGverts[iElement matches vertex name]
#           jIndex <- AGGverts[jElement matches vertext name]
#
#           W(jIndex,iIndex) <- 1/degree(i) (retrieved from degreeVect)
#
#       }
#
#============#Creating Matrix (I-(1-B)*W)============
#       I (identity matrix) <- matrix of dimensions (number of nodes x number of nodes), all zeroes
#
#       for (n in 1:sizeVerts (number of nodes))  {
#
#       I(n,n) <- 1
#
#       }
#
#       #Up to this point we have created the identity matrix
#
#
#       beta<-param[1]
#       matrix2 (matrix pre inversion) <- (I-(1-beta)*W)
#
# =========#Inverting matrix2=====================
#
#       inverted<-solve(matrix2)
#
#==========#Computing Matrix F=====================
#
#       F <- beta*inverted
#
#==========#Computing matrix E=====================
#
#       #create matrix h
#
#       h <- diag(Gene Score values, from AGGverts)
#
#       #post multiply F by h
#
#       E <- F*h
#
#==========#Extracting Scores from Matrix F========
#       influenceVect <- numeric(length = number of edges) #make a blank vector
#
#       AGGedges[,"Influence"]
#
#       for (n in 1:number of edges) {
#           iElement=AGGedges[n,from]
#           jElement=AGGedges[n,to]
#
#           #Find iElement and jElement in AGGverts to obtain their indices
#           iIndex <- AGGverts[iElement matches vertex name]
#           jIndex <- AGGverts[jElement matches vertext name]
#
#           Influence Score <- E(jIndex,iIndex)
#           AGGedges[n,"Influence"] <- Influence Score
#
#       #Up to here, we've assigned the influence score for an individual edge
#       }
#
#========Creating EGG==============
#
#
# EGGedges <- AGGedges
# EGGverts <- AGGverts
#
# EGG <- graph_from_data_frame(EGGedges,EGGverts)
#
#
#
# graphAttr(EGG,"EGGversion") <- EGGversion 1.0
#
#
#
# } #ending if statement that directs to Leiseron Method
#
# else {
#  stop("No method selected")
#}
#
# return(EGG)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
