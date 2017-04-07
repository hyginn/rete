#' Filter inconsistent CNA entries
#'
#' This filter removes genes from rCNA files that have inconsistent signs.
#' Raw rCNA files can contain genes whose copy numbers are consistently
#' increased or decreased relative to normal expression.  This filter finds
#' genes whose expression is not consistently increased or decreased, then
#' removes them from the dataset.
#'
#' @param rCNAfile The rCNA RDS file that contains the CNAs to consider
#' @param fNames A vector of filenames to filter
#' @param dOut The local path to a directory to output the filtered files
#' @param filterFile The filename to output the list of genes removed
#' @param iT The inconsistency threshold to use.  Default: 0.75
#' @param silent Whether or not to write progress information to console, default: FALSE
#' @param noLog Whether or not to log results, default: FALSE
#' @return The list of filtered genes
filter.inconsistentCNA <- function(rCNAfile, fNames, dOut, filterFile,
    iT=0.75, silent=FALSE, noLog=FALSE
) {
    # check for input file readability
    .filter.utils.fileReadable(rCNAfile)
    toFilter <- .filter.utils.fileReadable(fNames)

    # check for output directory writability
    .filter.utils.dirWritable(dOut)

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

    # Filter out selected genes 
    for (target in toFilter) {
        outputPath = .filter.utils.outputFilename(dOut, target)

        if (.filter.utils.getMagic(target) == 'MAF') {
            .filter.utils.filterrSNV(target, outputPath, genes)
        } else {
            .filter.utils.filterrCNA(target, outputPath, genes)
        }
    }

    return(genes)
}
