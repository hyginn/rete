# filter.utils.R

# Utility function for removing genes from a dataset
#
# Normalize the formatting of a tumor sample barcode.  Text files use
#   hyphens to separate elements, RDS converts that to '.'
#
# @param barcode A tumour barcode
# @param separator The separator used in this barcode, default '-'
# @return A normalized barcode that is upper case and has '.' as a separator
.filter.utils.normalizeBarcode <- function(
        barcode, separator="-") {

    if (typeof(barcode) != "character") {
        stop("Barcode must be of type character")
    }

    if (separator != ".") {
        barcode <- gsub('-', '.', barcode)
    }

    return(toupper(barcode))
}

# Extract the patient portion of a tumor sample barcode.  Expects a
#   normalized barcode.
#
# @param barcode A tumour barcode
# @param separator The separator used in this barcode, default '.'
# @return The barcode truncated to the study participant level,
#   with '.' as a separator
.filter.utils.patientFromBarcode <- function(
        barcode, separator=".") {

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

# Guess whether a file is MAF/rSNV or RDS/rCNA
#
# Gets the first byte of a file, return MAF if it starts with a hash
#   Not truly smart enough to handle all of our tab-delimited formats,
#   but it suffices for RDS vs (non-compressed) MAF.  If we were to 
#   be dealing with compressed MAF files, I could no longer recognize
#   the format based off the first byte of the file, and would instead
#   need to either wrap readRDS in a tryCatch() call or check the first
#   decompressed byte.  Decompressed RDS files do not start with a hash,
#   and MAF files do.
#
# @param filename The filename to check
# @return 'MAF' if the file starts with a hash, 'RDS' otherwise
.filter.utils.getMagic <- function(filename) {
    magic <- readBin(con=filename, what="int", n=1, size=1, signed=FALSE)

    if (magic == 35) {
        return("MAF")
    } else {
        return("RDS")
    }
}

# Generate a normalized output file path by concatenating the provided
#  output directory with the basename of the input file.
#
# @param path The desired output directory
# @param file A path whose base filename should be concatenated to the output
#   directory
# @return A canonical, concatenated path
.filter.utils.outputFilename <- function(path, file) {
    bName <- basename(file)
    outName <- normalizePath(paste(normalizePath(path), '/', bName, sep=''),
        mustWork=FALSE)
    return(outName)
}

# Validate basic file readability status
#
# @param fNames A vector of filenames to check for existence and readability
# @param level The level of error a non-readable file should hold.  If
#   error, execution stops.  Otherwise a warning to the logs is emitted.
# @return A vector of readable files
.filter.utils.fileReadable <- function(fNames, level='error') {
    # XXX: need tests
    readable <- c()
    for (fName in fNames) {
    	if (file.exists(fName) == FALSE) {
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
        readable <- c(readable, fName)
    }

    return(readable)
}

# Validate basic directory writability status
#
# @param dNames A vector of directory names to check for existence and
#   writeability
# @param level The level of error a non-writeable directory should hold.  If
#   error, execution stops.  Otherwise a warning to the logs is emitted.
# @return A vector of writeable directories
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
        writeable <- c(writeable, dName)
    }

    return(writeable)
}

# Scan an rSNV MAF file, outputting into the output file a new MAF
#   file that has only the records for which the filter says not to ignore
#
# @param inFile The input rSNV file
# @param outFile The output rSNV file
# @param removeGenes The list of genes to remove
# @return A two-element vector, where the first element is the number
#   of records kept from the original rSNV file and the second is
#   the number of records removed.
.filter.utils.filterrSNV <- function(inFile, outFile,
    removeGenes=c()) {

    rSNV <- file(inFile, open="r")
    SNV <- file(outFile, open="w")

    # Duplicate the header
    hasHeaderHash <- FALSE
    hasHeaderHugo <- FALSE
    line <- readLines(con=rSNV, n=1)
    while (length(line) != 0) {
        if (!startsWith(line[[1]], "#") &&
            !startsWith(line[[1]], "Hugo_Symbol")
        ) {
            break
        } else if (startsWith(line[[1]], '#')) {
            hasHeaderHash <- TRUE
        } else if (startsWith(line[[1]], 'Hugo_Symbol')) {
            hasHeaderHugo <- TRUE
        }
        writeLines(line, con=SNV)
        line <- readLines(con=rSNV, n=1)
    }
   
    if (!hasHeaderHash || !hasHeaderHugo) {
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
#
# @param inFile The input rCNA file
# @param outFile The output rCNA file
# @param removeGenes The list of genes to remove
# @return A two-element vector, where the first element is the number
#   of genes kept from the original rCNA file and the second is
#   the number of genes removed.
.filter.utils.filterrCNA <- function(inFile, outFile, removeGenes=c()) {
    rCNA <- readRDS(inFile)

    # Filter rCNA data
    removed <- 0
    keep <- c()
    for (gene in rownames(rCNA)) {
        if (gene %in% removeGenes) {
            removed <- removed + 1
            next
        }
        keep <- c(keep, gene)
    }
    CNA <- rCNA[keep, ]

    saveRDS(CNA, file=outFile)

    return(c(length(keep), removed))
}

# [END]
