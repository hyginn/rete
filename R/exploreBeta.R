# exploreBeta.R

#' Title.
#'
#' \code{exploreBeta} Function to find best beta parameter to use for heat diffusion 
#' of a graph.
#'
#' Details.
#' @section Mode: exploreBeta traverses the vertices and finds percentage of heat
#' that passes through the neighbouring vertices. The function travels a two edge
#' path and calculates the percentage of heat at that vertex. The function calculates 
#' the average heat at the neighbouring vertices and returns the smallest beta value
#' where no more than 5% of heat is passed along two edges away from the vertex where 
#' heat is applied. 
#'
#' @param <gG> <a graph of type igraph>.
#' 
#' @return <beta> a float number.
#'
#' @examples
#' exploreBeta(gG)
#' 
# Helper function to calculate heat
getF <- function(gG, B, v, F){
    # Find all neighbours 
    neighbours <- neighbors(gG, v) #Canadian eh!
    # Create adjecency matrix of subgraph of neighbours + vertex
    adjMatrix <- as_adj(induced_subgraph(gG, c(neighbours, v)))
    # Weight the matrix
    W <- adjMatrix/degree(gG, v)
    # Create identity matrix of equal degree 
    I <- (diag(1 + degree(gG, v)))
    # Return F
    return (F*(B * solve((I - ((1-B) * W)))))
    
}

explorebeta <- function(gG){
  
  # Paritiion Beta
  B <- c(0.05, 0.10, 0.15, 0.20, 0.25,
         0.30, 0.35, 0.40, 0.45, 0.50,
         0.55, 0.60, 0.65, 0.70, 0.75,
         0.80, 0.85, 0.95, 1.00)

  int = 1
  array = c()
  
  for(beta in B){
    nFluxAve = 0
    nFlux = 0
    gFluxAve = 0
    gFlux = 0
    
    for(v in 1:length(V(gG))){
      # Iterate over all vertices
      
      # Calculate heat (flux) assuming heat of degree 1 is applied to
      # vertex of interest
      flux <- getF(gG, beta, V(gG)[v], 1)
      value <- max(flux)
      # This "value" is the worst case heat that can be transferred
      neighbour <- neighbors(gG, V(gG)[v])
      
      for (n in 1:length(neighbour)){
        # Iterate over all neighbours of vertex to see how
        # much heat they will transfer to their neighbours
        flux <- getF(gG, beta, neighbour[n], value)
        # Worst case heat transfer is max of flux
        nFlux = nFlux + max(flux)
      }
      # Take average of heat over all neighbours
      nFluxAve = nFlux/length(neighbour) # Average the finals

      gFlux = gFlux + nFluxAve
    }
    # Take average of heat transwered by each vertex to 2 vertices away 
    gFluxAve = gFlux/length(V(gG))
    array[int] <- gFluxAve
    int = int + 1
  }
  # Return lowest possible beta that only allows for at most 5% heat transfer
  # Past 2 vertices away from vertex of interest throughout whole graph
  return(B[min(which(array < 0.05))])
}


                      