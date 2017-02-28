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

### Small utility functions.  Possibly of use elsewhere

# Extract the patient portion of a tumor sample barcode
.filter.unimportant_genes.patientFromBarcode <- function(
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

### Module-specific utility functions.

# Extract the barcode from the hash key
.filter.unimportant_genes.barcodeFromHashkey <- function(
        hashkey, separator="|") {

    if (typeof(hashkey) != "character") {
        stop("Hashkey must be of type character")
    }

    bits <- strsplit(hashkey, separator, fixed=TRUE)[[1]]
    if (length(bits) != 2) {
        stop("Hashkey incorrectly formatted")
    }

    if (length(bits[[2]]) == 0) {
        stop("Hashkey too short")
    }

    return(bits[[2]])
}

# Validate basic file status
.filter.unimportant_genes.validFiles <- function(
        exprData, clinData, fNames, dOut
    ) {
    # XXX: need tests
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

    toFilter <- c()
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

        toFilter <- c(toFilter, filename)
    }

    return(toFilter)
}

# Load the expression rCNA file
.filter.unimportant_genes.loadExpression <- function(filename, rT) {
    genehash <- new.env(hash=TRUE)
    numSampleReads <- 0

    samples <- readRDS(filename)
    for (gene in rownames(samples)) {
        for (s in colnames(samples)) {
            # Preset everything for the sample's hashtable record
            df <- list(ignore = FALSE, numReads = 0)

            df$numReads <- df$numReads + samples[gene, s]

            if (df$numReads > rT) {
                numSampleReads <- numSampleReads + 1
            }

            # Store the number of reads
            genehash[[paste(gene, s, sep='|')]] <- df
        }
    }

    return(list("numSampleReads"=numSampleReads, "samples"=genehash))
}

# Filter the expression RDS data
.filter.unimportant_genes.filterExpression <- function(
        filter, survivalTime, rT, pT, cT
    ) {

    samples <- filter$samples
    numSampleReads <- filter$numSampleReads

    # Operate over all elements in the environment hash filter
    #   This is /not/ a copy.  Changes to filter persist are
    #   visible out of this scope.
    for (s in ls(samples)) {
        numReads <- samples[[s]]$numReads

        if (numReads < rT) {
            samples[[s]]$ignore <- TRUE
        } else if (numSampleReads < (pT * numReads)) {
# XXX: Still need test for this case
            samples[[s]]$ignore <- TRUE
        }

        # Note that survivalTime keys are not the full length
        #   of a full sample barcode key.
        patient <- toupper(.filter.unimportant_genes.patientFromBarcode(
            .filter.unimportant_genes.barcodeFromHashkey(s),
            separator="."))
        if (survivalTime[[patient]] < cT) {
            samples[[s]]$ignore <- FALSE
        }
    }

    return(samples)
}

# Load the clinical data RDS file, text format?
.filter.unimportant_genes.loadVitalStatus <- function(filename) {
    # XXX: need tests
    readhash <- new.env(hash=TRUE)

    clin <- readRDS(filename) # XXX: is this really an RDS file?

    # Select cases to consider
    patients <- colnames(clin)
    for (patient in patients) {
        vs <- clin["vital_status", patient]
        survivalTime <- -1

        if (vs == 1) { # dead
            survivalTime <- clin["days_to_death", patient]
        } else if (vs == 0) { # alive
            dtd <- clin["days_to_death", patient]
            if (dtd != NA) {
                next # exclude, a living dead person
            }

            dtlf <- clin["days_to_last_followup", patient]
            if (dtlf != NA) {
                survivalTime <- dtlf
            }
        } else if (vs == NA) { # don't know
            next # exclude
        } else {
            stop("Invalid vital_status value in clinical data")
        }

        readhash[[patient]] <- survivalTime
    }

    return(readhash)
}

# Scan an rMUT MAF file, outputting into the output directory a new MAF
#   file that has only the records for which the filter says not to ignore
.filter.unimportant_genes.processrMUT <- function(filename, dOut, filter) {
    outName <- .filter.unimportant_genes.outputFilename(dOut, filename)

    rMUT <- file(filename, open="r")
    MUT <- file(outName, open="w")

    # Duplicate the header
    hasHeader <- FALSE
    line <- readLines(con=rMUT, n=1)
    while (length(line) != 0) {
        if (!startsWith(line[[1]], "#")) {
            break
        }
        hasHeader <- TRUE
        writeLines(line, con=MUT)
        line <- readLines(con=rMUT, n=1)
    }
   
    if (!hasHeader) {
        stop("Malformed rMUT file: missing header")
    }
 
    # Filter the data
    while (length(line) != 0) {
        fields <- strsplit(line[[1]], "\t", fixed=TRUE)[[1]]
        if (length(fields) != 12) {
            stop("Malformed rMUT file: incorrect field count")
        }

        record <- paste(c(fields[1],
            gsub('-', '.', fields[11])), collapse="|")
        if (exists(record, envir=filter) &&
            filter[[record]]$ignore == TRUE) {
            line <- readLines(con=rMUT, n=1)
            next
        }
        writeLines(line, con=MUT)
        line <- readLines(con=rMUT, n=1)
    }

    close(MUT)
    close(rMUT)

    return(outName)
}

# Scan an rCNA RDS object, outputting into the output directory a new RDS
#   file that has only the records for which the filter says not to ignore
.filter.unimportant_genes.processrCNA <- function(filename, dOut, filter) {
    # XXX: need tests
    outName <- .filter.unimportant_genes.outputFilename(dOut, filename)

    rCNA <- readRDS(filename)

    # Filter rCNA data
    samples <- intersect(colnames(rCNA), ls(filter))
    keep <- c()
    for (gene in genes) {
        if (filter[[gene]]$ignore == FALSE) {
             keep <- c(keep, gene)
        }
    }
    CNA <- rCNA[, keep]

    saveRDS(rCNA, file=outName)
}

### Main function
filter.unimportant_genes <- function(exprData, clinData, fNames, dOut,
        rT=3, pT=0.7, cT=5, silent = FALSE, noLog = FALSE
    ) {
    # XXX: need tests
    # Validate file existence
    toFilter <- .filter.unimportant_genes.validFiles(
        exprData, clinData, fNames, dOut)

    # Count the reads per gene and sample in the gene expression file
    expression <- .filter.unimportant_genes.loadExpression(exprData, rT)

    # Survival time data
    survivalTime <- .filter.unimportant_genes.loadVitalStatus(clinData)
    
    # Filter the expression data
    expression <- .filter.unimportant_genes.filterExpression(
        readhash, expression, survivalTime, rT, pT, cT)
    # expression is now an environment hash

    # Process rMUT and rCNA files
    for (filename in toFilter) {
        fType <- .filter.unimportant_genes.getMagic(filename)
        if (fType  == 'MAF') {
            rMUT <- .filter.unimportant_genes.processrMUT(
                filename, dOut, expression)
        } else {
            rCNA <- .filter.unimportant_genes.processrCNA(
                filename, dOut, expression)
        }
    }

    return(NA)
}

# [END]
