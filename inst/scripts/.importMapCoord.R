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

makeMap <- function(ucscMod) {
    if(length(ucscMod$exonStart) != length(ucscMod$exonEnd)) {
        return(NULL)
    }

    # prepare the environment (exon matrix, start, end)
    exons <- cbind(ucscMod$exonStart, ucscMod$exonEnd)
    start <- ucscMod$start
    end <- ucscMod$end
    if(is.unsorted(c(start, as.vector(t(exons)), end))) {
        return(NULL)
    }

    # return a function that uses the exon matrix and start, stop positions
    return(
        function(x) {
            if (x < start | x > end) {
                # out-of-range coordinate
                return(NULL)
            }
            iEx <- which(x >= exons[ , 1] & x <= exons[ , 2])
            if(length(iEx) != 1){
                # intron
                return(1)
            }

            if (iEx == 1) {
                lPre <- 0
            } else {
                lPre <- sum((exons[ , 2] - exons [ , 1] + 1)[1:(iEx-1)])
            }
            lIn <- x - exons[iEx, 1] + 1

            return(lPre + lIn)
        }
    )
}

mapCoord <- function(geneData, dropExonBoundaries = TRUE, silent = FALSE, writeLog = TRUE){
    library(parallel)
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores, type = "PSOCK")

    # "function factory" manufacturing step, parallelized
    geneData$map <- parLapply(geneData, function(x) makeMap(x))

    saveRDS(geneData, file = "genedata.rds")

    stopCluster(cl)
}



# [END]


# Validation of makeMap (function for one gene model )

allSamples <- list()
allSamples$Bob <- list(name="bob", strand='+', start=200, end=300, exonStart=c(201,220,240), exonEnd=c(210,230,250))

sampleFx <- makeMap(allSamples$Bob)

sampleFx(202)
sampleFx(222)
sampleFx(242)



