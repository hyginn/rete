# filter.UNIMPORTANT_GENES.R

#' Remove "unimportant genes" from a dataset
#'
#' \code{filter.unimportant_genes} placeholder function
#'
#' @param exprData The local fully qualified filename of the gene expression data RDS file.
#' @param clinData The local fully qualified filename of the clinical data RDS file.
#' @param fNames A vector of local, fully-qualified filenames of rMUT and/or rCNA RDS files.
#' @param dOut The local, fully-qualified path of a directory to store the output rMUT and/or rCNA RDS files.
#' @param rT Threshold: only samples with numReads < rT are retained in the output.
#' @param pT Threshold: only samples with numSampleReads < (pT * number of samples) are retained in the output.
#' @param cT Threshold: only samples with "days_to_death_or_fup" < cT are retained in the output.
#' @param silent Boolean. Default: FALSE. Whether or not to write progress information to the console.
#' @param noLog Boolean. Default: FALSE. Whether or not to log results.
#' @export

# Load the expression RDS file
.filter.unimportant_genes.loadExpression <- function(filename, vital_status, rT, pT, cT) {
    readhash <- new.env(hash=TRUE)

    stored <- readRDS(filename)
    for (foo in stored) {
        df = data.frame(ignore = FALSE, numReads = foo$reads)

        readhash[[foo$gene]] <- df
    }

    return readhash
}

# Load the clinical data RDS file
.filter.unimportant_genes.loadVitalStatus <- function(filename) {
    readhash <- new.env(hash=TRUE)

    stored <- readRDS(filename)

    return readhash
}

# Get the first byte of a file, return MAF if it starts with a hash
#   Not truly smart enough to handle all of our tab-delimited formats,
#   but it suffices for RDS vs (non-compressed) MAF.  If we were to 
#   be dealing with compressed MAF files, I could no longer recognize
#   the format based off the first byte of the file, and would instead
#   need to either wrap readRDS in a tryCatch() call or check the first
#   decompressed byte.  Decompressed RDS files do not start with a hash,
#   and MAF files do.
.filter.unimportant_genes.getMagic <- function(filename) {
    magic <- readBin(con=filename, what="int", n=1, size=1, signed=FALSE)

    if (magic == 35) {
        return("MAF")
    } else {
        return("RDS")
    }
}

# Generate a normalized output file path by concatenating the provided
#   output directory with the basename of the input file.
.filter.unimportant_genes.outputFilename <- function(path, file) {
    bName <- basename(file)
    outName <- normalizePath(paste(normalizePath(path), '/', bName, sep=''),
        mustWork=FALSE)
    return(outName)
}

# Scan an rMUT MAF file, outputting into the output directory a new MAF
#   file that has only the records for which the filter says not to ignore
.filter.unimportant_genes.processrMUT <- function(filename, dOut, filter) {
    outName <- .filter.unimportant_genes.outputFilename(dOut, filename)

    rMUT <- file(filename, "r")
    close(rMUT)
}

# Scan an rCNA RDS object, outputting into the output directory a new RDS
#   file that has only the records for which the filter says not to ignore
.filter.unimportant_genes.processrCNA <- function(filename, dOut, expression) {
    outName <- .filter.unimportant_genes.outputFilename(dOut, filename)

    rCNA <- readRDS(filename)
    saveRDS(rCNA, file=outName)
}

# Main function
filter.unimportant_genes <- function(exprData, clinData, fNames, dOut,
        rT, pT, cT, silent = FALSE, noLog = FALSE
    ) {
    # Validate file existence
    if (file.exists(exprData) == FALSE) {
        stop(paste('ERROR: Expression data file', exprData, 'does not exist.'))
    }
    if (file.access(exprData, mode=4) == -1) {
        stop(paste('ERROR: Expression data file', exprData, 'is not readable.'))
    }
    if (file.exists(clinData) == FALSE) {
        stop(paste('ERROR: Clinical data file', clinData, 'does not exist.'))
    }
    if (file.access(clinData, mode=4) == -1) {
        stop(paste('ERROR: Clinical data file', clinData, 'is not readable.'))
    }
    if (file.exists(dOut) == FALSE) {
        stop(paste('ERROR: Output directory', dOut, 'does not exist.'))
    }
    if (file.info(dOut)$isdir == FALSE) {
        stop(paste('ERROR: Output directory', dOut, 'is not a directory.'))
    }
    if (file.access(dOut, mode=2) == -1) {
        stop(paste('ERROR: Output directory', dOut, 'is not writable.'))
    }

    # Survival time data
    vitalStatus <- .filter.unimportant_genes.loadVitalStatus(clinData)

    # Count the reads per gene and sample in the gene expression file
    # while filtering for expression
    expression <- .filter.unimportant_genes.loadExpression(exprData, vitalStatus, rT, pT, cT)

    # Process rMUT and rCNA files
    for (filename in fNames) {
        if (file.exists(filename) == FALSE) {
            .appendToLog(paste('WARN: Data file', filename,
                'does not exist, skipping.'))
            next
        }
        if (file.access(filename, mode=4) == -1) {
            .appendToLog(paste('WARN: Data file', filename,
                'is not readable, skipping.'))
            next
        }

        fType <- .filter.unimportant_genes.getMagic(filename)
        if (fType  == 'MAF') {
            rMUT <- .filter.unimportant_genes.processrMUT(filename, dOut, expression)
        } else {
            rCNA <- .filter.unimportant_genes.processrCNA(filename, dOut, expression)
        }
    }

    return(TRUE)
}

# [END]
