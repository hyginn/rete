# <non-exported function header>
# <check best practices>

mapCoord <- function(geneData, dropExonBoundaries = TRUE, silent = FALSE, writeLog = TRUE){
  library(parallel)
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores, type = "PSOCK")
    
  
  #might not work 
  geneData$map <- parLapply(geneData, 1, function(x) makeMap(x))
  
  saveRDS(geneData, file = "example_data_with_functions.RDS")
  
  stopCluster(cl)
}

makeMap <- function(gMod) {
  # prepare the exon matrix
  lG <- length(gMod)
  exons <- matrix(numeric(lG - 2), ncol = 2)
  exons[, 1] <- (gMod[2:(lG-1)])[c(TRUE, FALSE)]
  exons[, 2] <- (gMod[3:(lG)])[c(TRUE, FALSE)]

  # return a function that uses the exon matrix
  return(
    function(x) {
      iEx <- which(x >= exons[ , 1] & x <= exons[ , 2])
      if(length(iEx) == 0){ return(0) }
      lPre <- sum((exons[ , 2] - exons [ , 1])[1:(iEx-1)])
      lIn <- x - exons[iEx, 1] + 1
      return(lPre + lIn)
    }
  )
}



