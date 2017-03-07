#' Filter inconsistent CNA entries
#'
#' This filter removes genes from rCNA files that have inconsistent signs.
#' Raw rCNA files can contain genes whose copy numbers are consistently
#' increased or decreased relative to normal expression.  This filter finds
#' genes whose expression is not consistently increased or decreased, then
#' removes them from the dataset.
#'
#' @param rCNA An rCNA RDS file to operate over
#' @param CNAFile The filename to output the filtered CNA data into
#' @param filterFile The filename to output the list of genes removed
#' @param iT The inconsistency threshold to use.  Default: 0.75
#' @return The list of filtered genes
filter.inconsistentCNA <- function(rCNAfile, CNAFile, filterFile,
    iT=0.75
) {
    # check for input file readability
    .filter.utils.fileReadable(rCNAfile)

    # load input data
    rCNA <- readRDS(rCNAfile)

    # generate filter list
    genes <- c()
    samples <- colnames(rCNA)
    sampleCount <- ncol(rCNA)
    for (gene in rownames(rCNA)) {
        amp <- 0
        del <- 0
        for (sample in samples) {
            val <- rCNA[gene, sample]
            if (val < 0) {
                del <- del + 1
            } else if (val > 0) {
                amp <- amp + 1
            }
        }

        if (del/sampleCount < iT && amp/sampleCount < iT) {
            genes <- c(genes, gene)
        }
    }

    # Save filter list
    writeLines(genes, con=filterFile)

    # Filter rCNA data
    .filter.utils.filterrCNA(rCNAfile, CNAFile, genes)
    return(genes)
}
