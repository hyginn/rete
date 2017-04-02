# <non-exported function header>

# mapCoord.R

#' Title.
#'
#' \code{<function>} description.
#'
#' Details.
#' @section <title>: Additional explanation.
#'
#' @param <p> <description>.
#' @param <q> <description>.
#' @return <description>.
#'
#' @family <optional description of family>
#'
#' @seealso \code{\link{<function>}} <describe related function>, ... .
#'
#' @examples
#' multiply(1, 1)
#' multiply(1i, 1i)
#' multiply(13, 3)
mapCoord <- function(geneData, dropExonBoundaries = TRUE, silent = FALSE, writeLog = TRUE){
    library(parallel)
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores, type = "PSOCK")

    makeMap <- function(ucscMod) {

        # check for valid data
        if(length(ucscMod$ExonStart) != length(ucscMod$ExonEnd)) { return(0) }

        # prepare the exon matrix
        exons <- cbind(ucscMod$ExonStart, ucscMod$ExonEnd)

        # return a function that uses the exon matrix
        return(
            function(x) {
                iEx <- which(x >= exons[ , 1] & x <= exons[ , 2])
                if(length(iEx) != 1){ return(0) }
                lPre <- sum((exons[ , 2] - exons [ , 1])[1:(iEx-1)])
                lIn <- x - exons[iEx, 1] + 1
                return(lPre + lIn)
            }
        )
    }

    # "function factory" manufacturing step, parallelized
    geneData$map <- parLapply(geneData, function(x) makeMap(x))

    saveRDS(geneData, file = "genedata.rds")

    stopCluster(cl)
}

# [END]


# Validation of makeMap (function for one gene model )


sample <- list(name="bob", start = 200, end = 300)
sample$ExonStart <- c(201,220,240)
sample$ExonEnd <- c(210,230,250)



