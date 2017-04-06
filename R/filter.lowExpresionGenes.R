#' Filter low expression genes
#'
#' This filter removes genes from rCNA and rSNV files that are not highly
#' expressed in the provided RNAseq files.
#'
#' @param exprData The normalized Firehose RNAseq files to scan for expression levels
#' @param rCNA The rCNA RDS files to filter out lower expression genes
#' @param rCNA The rSNV text files to filter out lower expression genes
#' @param dOut The local path to a directory to output the filtered files
#' @param whitelist The name of the file that contains genes to preserve
#' @param rT The minimum number of reads for genes to keep.  Default: 3
#' @param pT The percentage of samples that must have at least rT reads.  Default: 0.7
#' @param silent Whether or not to write progress information to console, default: FALSE
#' @param noLog Whether or not to log results, default: FALSE
#' @return The list of filtered genes
filter.lowExpressionGenes <- function(exprData, rCNA, rSNV, dOut, whitelist,
    rT=3, pT=0.7,
    silent=FALSE, noLog=FALSE
) {
    # check numeric arguments
    err <- .checkArgs(rT, like = numeric(1))
    if (length(err) > 0) {
        stop(err)
    }
    err <- .checkArgs(pT, numeric(1))
    if (length(err) > 0) {
        stop(err)
    } else if (pT > 1 || pT < 0) {
        stop('pT must be in the range of 0.0 - 1.0')
    }

    # check file arguments
    err <- .checkArgs(dOut, 'DIR')
    if (length(err) > 0) {
        stop(err)
    }

    for (fName in c(rCNA, rSNV, whitelist, exprData)) {
        err <- .checkArgs(fName, 'FILE_E')
        if (length(err) > 0) {
            stop(err)
        }
    }

    # check for input file readability
    .filter.utils.fileReadable(exprData)
    .filter.utils.fileReadable(whitelist)
    .filter.utils.fileReadable(rSNV)
    .filter.utils.fileReadable(rCNA)

    # generate pre-whitelist list of genes to remove
    rawRemove <- .filter.lowExpressionGenes.processExpressionData(
        exprData, rT, pT)

    # read whitelist file
    toKeep <- readLines(con=whitelist)

    # remove whitelist entries from removal list
    toRemove <- rawRemove[!(rawRemove %in% toKeep)]

##    # save filter list
##    filter <- tempfile()
##    writeLines(toRemove, filter)

    currDir <- getwd()
    setwd(dOut)

    # Filter out selected genes 
    for (target in rCNA) {
        fName <- basename(target)
        .filter.utils.filterrCNA(target, fName, toRemove)
    }
    for (target in rSNV) {
        fName <- basename(target)
        .filter.utils.filterrSNV(target, fName, toRemove)
    }

    # cleanup files
    setwd(currDir)
##    unlink(filter)

    return(toRemove)
}

# Parse the expression data
.filter.lowExpressionGenes.parseExpressionFile <- function(filename) {
    expr <- file(filename, open="r")
   
    # Get the header
    line <- readLines(con=expr, n=2)
    if (length(line) == 0) {
        close(expr)
        stop("Expression data file empty")
    } else if (length(line) == 1) {
        close(expr)
        stop("Expression data headers incomplete")
    }

    # Save the barcode ids
    barcodes <- strsplit(line, "\t")[[1]]
    colcount <- length(barcodes)

    # Check for data beyond gene names
    if (colcount < 2) {
        close(expr)
        stop("Expression data missing samples")
    }

    # Only consider tumour samples
    keepList <- .filter.utils.sampleTypeFromBarcode(
        .filter.utils.normalizeBarcode(barcodes[2:colcount]))$type == 'T'

    # Check the data type
    dataType <- strsplit(line, "\t")[[2]][2]
    if (dataType != 'normalized_count') {
        close(expr)
        stop("Expression data type incorrect")
    }
   
    data <- data.frame(matrix(nrow=0, ncol=sum(keepList)),
        stringsAsFactors=FALSE)
    line <- readLines(con=expr, n=1)[1]
    while (length(line) != 0) {
        fields <- strsplit(line, "\t")[[1]]
        gene <- strsplit(fields[1], '|', fixed=TRUE)[[1]][1]

        samples <- suppressWarnings(as.numeric(fields[2:colcount]))
        if (any(is.na(samples))) {
            close(expr)
            stop("Expression data non-numeric")
        }

        data[gene, ] <- samples[keepList]

        line <- readLines(con=expr, n=1)
    }

    colnames(data) <- make.names(toupper(barcodes[c(FALSE, keepList)]))

    close(expr)
    return(data)
}

# Load the expression data, return list of low-expression genes
.filter.lowExpressionGenes.processExpressionData <- function(files, rT, pT) {
    # data has not been log2 transformed, untransform threshold
    exp2 <- 2^rT

    summary <- data.frame(matrix(nrow=0, ncol=2))
    colnames(summary) <- c('expressed', 'samples')

    # scan the RNAseq files
    for (filename in files) {
        samples <- .filter.lowExpressionGenes.parseExpressionFile(filename)
        samplesInFile <- length(colnames(samples))
        for (gene in rownames(samples)) {
            expressed <- 0

            for (barcode in colnames(samples)) {
                # Count the number of samples that express this gene
                if (samples[gene, barcode] > exp2) {
                    expressed <- expressed + 1
                }
            }

            # Merge this data in with the existing summmary
            oldSummary <- summary[gene,]
            if (is.na(oldSummary$expressed)) {
                summary[gene, ] <- c(expressed, samplesInFile)
            } else {
                summary[gene, ] <- c(
                    oldSummary$expressed+expressed,
                    oldSummary$samples+samplesInFile)
            }
        }
    }

    # calculate the percentage expressed
    lowExpression <- rownames(summary[summary$expressed/summary$samples < pT,])

    return(lowExpression)
}
