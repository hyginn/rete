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

PSEUDOCODE

exploreBeta <- function(gG){
possibleBeta <- c(0-1 indices)
for beta in possiblebeta
    for vertex in gG
      calculate heat transfusion on vertex
      find neighbours
      for neighbours of vertex
          find second_neighbours
          calculate heat transfusion to second_neighbours
      calculate average heat from vertex to second neihgbours
    calculate average heat of all vertices to their second neighbours
  return lowest beta with average heat from vertex to second neighbours loess than 5%

}

# [END]