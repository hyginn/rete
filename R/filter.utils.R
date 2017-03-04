# filter.utils.R

#' Remove "unimportant genes" from a dataset
#'
#' \code{filter.utils} placeholder function
#'
#' @param exprData The local fully qualified filename of the gene expression data RDS file.
#' @param clinData The local fully qualified filename of the clinical data RDS file.
#' @param fNames A vector of local, fully-qualified filenames of rSNV and/or rCNA RDS files.
#' @param dOut The local, fully-qualified path of a directory to store the output rSNV and/or rCNA RDS files.
#' @param rT Threshold: only samples with numReads < rT are retained in the output.
#' @param pT Threshold: only samples with numSampleReads < (pT * number of samples) are retained in the output.
#' @param cT Threshold: only samples with "days_to_death_or_fup" < cT are retained in the output.
#' @param silent Boolean. Default: FALSE. Whether or not to write progress information to the console.
#' @param noLog Boolean. Default: FALSE. Whether or not to log results.

### Small utility functions.  Possibly of use elsewhere

# Extract the patient portion of a tumor sample barcode
.filter.utils.patientFromBarcode <- function(
        barcode, separator="-") {

    if (typeof(barcode) != "character") {
        stop("Barcode must be of type character")
    }

    bits <- strsplit(barcode, separator, fixed=TRUE)[[1]]
    if (length(bits) < 3) {
        stop("Barcode does not contain a patient identifier")
    }

    patient <- paste(bits[1:3], collapse=separator)

    return(patient)
}

# Get the first byte of a file, return MAF if it starts with a hash
#   Not truly smart enough to handle all of our tab-delimited formats,
#   but it suffices for RDS vs (non-compressed) MAF.  If we were to 
#   be dealing with compressed MAF files, I could no longer recognize
#   the format based off the first byte of the file, and would instead
#   need to either wrap readRDS in a tryCatch() call or check the first
#   decompressed byte.  Decompressed RDS files do not start with a hash,
#   and MAF files do.
.filter.utils.getMagic <- function(filename) {
    magic <- readBin(con=filename, what="int", n=1, size=1, signed=FALSE)

    if (magic == 35) {
        return("MAF")
    } else {
        return("RDS")
    }
}

# Generate a normalized output file path by concatenating the provided
#   output directory with the basename of the input file.
.filter.utils.outputFilename <- function(path, file) {
    bName <- basename(file)
    outName <- normalizePath(paste(normalizePath(path), '/', bName, sep=''),
        mustWork=FALSE)
    return(outName)
}

# Validate basic file readability status
.filter.utils.fileReadable <- function(fNames, level='error') {
    # XXX: need tests
    readable <- c()
    for (fName in fNames) {
    	if (file.exists(exprData) == FALSE) {
            if (level == 'error') {
    	        stop(paste('ERROR: filter data file', fName, 'does not exist.'))
    	    } else {
                .appendToLog(paste('WARN: filter data file', fName,
                    'does not exist, skipping.'))
                next
    	    }
    	}
    	if (file.access(fName, mode=4) == -1) {
            if (level == 'error') {
    	        stop(paste('ERROR: filter data file', fName,
                    'is not readable.'))
    	    } else {
                .appendToLog(paste('WARN: filter data file', fName,
                    'is not readable, skipping.'))
                next
    	    }
    	}
        readable <- c(readable, filename)
    }

    return(readable)
}

# Validate basic directory writability status
.filter.utils.dirWritable <- function(dNames, level='error') {
    # XXX: need tests
    writeable <- c()
    for (dName in dNames) {
    	if (file.exists(dName) == FALSE) {
            if (level == 'error') {
    	        stop(paste('ERROR: directory', dName, 'does not exist.'))
    	    } else {
                .appendToLog(paste('WARN: directory', dName,
                    'does not exist, skipping.'))
                next
    	    }
    	}
    	if (file.info(dName)$isdir == FALSE) {
            if (level == 'error') {
    	        stop(paste('ERROR: directory', dName, 'is not a directory.'))
    	    } else {
                .appendToLog(paste('WARN: directory', dName,
                    'is not a directory, skipping.'))
                next
    	    }
    	}
    	if (file.access(dName, mode=2) == -1) {
            if (level == 'error') {
    	        stop(paste('ERROR: directory', dName,
                    'is not writeable'))
    	    } else {
                .appendToLog(paste('WARN: directory', dName,
                    'writeable, skipping.'))
                next
    	    }
    	}
        writeable <- c(writeable, filename)
    }

    return(writeable)
}

# Scan an rSNV MAF file, outputting into the output directory a new MAF
#   file that has only the records for which the filter says not to ignore
filter.utils.processrSNV <- function(inFile, outFile,
    removeGenes=c()) {

    rSNV <- file(inFile, open="r")
    SNV <- file(outFile, open="w")

    # Duplicate the header
    hasHeader <- FALSE
    line <- readLines(con=rSNV, n=1)
    while (length(line) != 0) {
        if (!startsWith(line[[1]], "#")) {
            break
        }
        hasHeader <- TRUE
        writeLines(line, con=SNV)
        line <- readLines(con=rSNV, n=1)
    }
   
    if (!hasHeader) {
        stop("Malformed rSNV file: missing header")
    }
 
    # Filter the data
    removed <- 0
    kept <- 0
    while (length(line) != 0) {
        fields <- strsplit(line[[1]], "\t", fixed=TRUE)[[1]]
        if (length(fields) != 12) {
            stop("Malformed rSNV file: incorrect field count")
        }

        # RDS converts the "-" in the barcodes to "."
        #    fields[11] %in% removeSamples) {
        if (fields[1] %in% removeGenes) {
            removed <- removed + 1
            line <- readLines(con=rSNV, n=1)
            next
        }
        kept <- kept + 1
        writeLines(line, con=SNV)
        line <- readLines(con=rSNV, n=1)
    }

    close(SNV)
    close(rSNV)

    return(c(kept, removed))
}

# Scan an rCNA RDS object, outputting into the output directory a new RDS
#   file that has only the records for which the filter says not to ignore
filter.utils.processrCNA <- function(inFile, outFile, removeGenes=c()) {
    rCNA <- readRDS(inFile)

    # Filter rCNA data
    removed <- 0
    keep <- c()
    for (gene in rownames(rCNA)) {
        if (gene %in% removeGenes) {
            removed <- removed + 1
            break
        }
        keep <- c(keep, gene)
    }
    CNA <- rCNA[, keep]

    saveRDS(CNA, file=outFile)

    return(c(length(keep), removed))
}

# [END]
