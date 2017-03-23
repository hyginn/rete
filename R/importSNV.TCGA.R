# importSNV.TCGA.R

#' Import MAF files from source and converts them to rSNV files.
#'
#' \code{importSNV.TCGA} Imports data from one or more files of MAF format, creates rSNV or adds to
#' existing rSNV file, returns rSNV file.
#'
#' @section Validations: validations sections has 3 functions - function that checks MAF file extension,
#' function that checks MAF file header, function that checks rSNV file header.
#'
#' @section Read MAF file: maf file is read by extracting require columns from maf file followed
#' by checking if every item in the class, type and chr columns are from a known and required list of values.
#' These columns as a data table are then added to an rSNV file.
#'
#' @section Wrtie Log: log section writes file processing progress along with other errors
#' that may occur.
#'
#' @param fMAF vector of MAF file names.
#' @param rSNV rSNV file to merged new MAF data with.
#' @param na.rm Remove all rows with missing or "NA" values and drop them from data frame.
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#'
#' @family
#'
#'   ## @seealso \code{\link{importM}} importM() is used to import ...
#'
#'   ## @examples ## \dontrun {
#'   ## importSNV.TCGA(c("aMAFfile.maf","bMAFfile"), "sampleSNV", na.rm = FALSE) ## }
#' @export

importSNV.TCGA <- function(  fMAF,
                             rSNV,
                             na.rm = FALSE,
                             silent = FALSE,
                             writeLog = TRUE) {

    myNotes <- character()
    myCall <- character()

    # ==== PARAMETERS ==========================================================

    # NL <- .PlatformLineBreak()
    # # Read header and one contents line
    # tmp <- readLines(fName, n = 2)

    # ==== VALIDATIONS =========================================================

    # General parameter checks
    cR <- character()
    cR <- c(cR, .checkArgs(fMAF,         like = rep("FILE_E", length(fMAF)),
                                                             checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,     like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(na.rm,        like = logical(1), checkSize = TRUE))

    if(length(cR) > 0) {
        stop(cR)
    }

    # ==== Validata rSNV  =======================================
    if (missing(rSNV)) {
        rSNVu <- "MAFtoSNV"

        write("sym\tchr\tstart\tend\tstrand\tclass\ttype\taRef\ta1\ta2\tTCGA\tUUID",
              rSNVu, append=FALSE)

        snvValidity <- character()
    } else {
        rSNVu <- rSNV
        rSNVf <- read.table(rSNV, header = TRUE)[,1]
        rowsBefore <- length(rSNVf)
        myNotes <- c(myNotes, paste("Number of rows originally in rSNV - ", rowsBefore))
        snvValidity <- .rSNVheadercheck(rSNV)
    }

    # ==== READ DATA ===========================================================
    filesFinished <- 0
    numFilesInvalid <- 0

    for (mafFile in fMAF) {
        if (length(.fileValidity(mafFile)) + length(.mafHeaderCheck(mafFile)) + length(snvValidity) == 0) {

            fileasDT <- read.table(mafFile, header = TRUE, fill = NA,  stringsAsFactors = FALSE,
                                        na.strings = " ")[,c(1,5,6,7,8,9,10,11,12,13,16,33)]

            # check if variant_class contains valid values from known list

            cat(sprintf("Processing file - \"%s\" - (%i/%i)", mafFile, filesFinished, length(fMAF)))

            classValid <- TRUE
            typeValid <- TRUE
            chrValid <- TRUE

            classVknown <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                        "Nonsense_Mutation", "Silent", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation",
                        "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR1" , "Intron", "RNA", "Targeted_Region")
            myVclasses <- unique(fileasDT[,6])

            if (! all(myVclasses %in% classVknown)) {
                # drop file
                classValid <- FALSE
                myNotes <- c(myNotes,paste(mafFile, " - was skipped because of unrecognized string in variant class column."))
            }

            # check if vairant_type contains valid values from known list
            typeVknown <- c("SNP", "DNP", "TNP", "ONP", "INS", "DEL", "Consolidated2")
            myVtypes <- unique(fileasDT[,7])
            if (! all(myVtypes %in% typeVknown)) {
                # drop file
                typeValid <- FALSE
                myNotes <- c(myNotes, paste(mafFile, " - was skipped because of unrecognized string in variant type column."))
            }

            # check chromosome value is valid
            typeChr <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
                         "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")
            myVchr <- unique(fileasDT[,2])
            if (! all(myVchr %in% typeChr)) {
                # drop file
                chrValid <- FALSE
                myNotes <- c(myNotes, paste(mafFile, " - was skipped because of unrecognized string in variant chromosome column."))
            }

            if (! classValid && typeValid && chrValid) {
                # there was an invalid value
                # drop file, go to next
                numFilesInvalid <- numFilesInvalid + 1
            } else {
                # all tests passed
                # add data to rSNV
                if (na.rm) {
                    fileasDT <- na.omit(fileasDT[,6])
                }

                write.table(fileasDT,  file=rSNVu, sep="\t", append=TRUE, row.names =
                                FALSE, col.names = FALSE, quote = FALSE)
            }
            filesFinished <- filesFinished + 1
        } else {
            numFilesInvalid <- numFilesInvalid + 1
            myNotes <- c(myNotes, .fileValidity(mafFile))
            myNotes <- c(myNotes, .mafHeaderCheck(mafFile))
            myNotes <- c(myNotes, snvValidity)
        }


        # ==== WRITE LOG ===========================================================

        if(writeLog) {

            myTitle <- "importSNV.TCGA"

            # Compile function call record
            if (missing(rSNV)) {
                myCall[1] <- "importSNV.TCGA("
                myCall[2] <- sprintf("fMAF = \"%s\", ", fMAF)
                myCall[3] <- sprintf("silent = %s, ", as.character(silent))
                myCall[4] <- sprintf("writeLog = %s, ", as.character(writeLog))
                myCall[5] <- sprintf("na.rm = %s)", as.character(na.rm))
                myCall <- paste0(myCall, collapse = "")
            } else {
                myCall[1] <- "importSNV.TCGA("
                myCall[2] <- sprintf("fMAF = \"%s\", ", fMAF)
                myCall[3] <- sprintf("rSNV = \"%s\", ", rSNVu)
                myCall[4] <- sprintf("silent = %s, ", as.character(silent))
                myCall[5] <- sprintf("writeLog = %s, ", as.character(writeLog))
                myCall[6] <- sprintf("na.rm = %s)", as.character(na.rm))
                myCall <- paste0(myCall, collapse = "")
            }

            # indicate input object name(s)
            #   (NA for this importSNV.TCGA())

            # Record progress information
            myNotes <- c(myNotes, paste("Total number of files in input fMAF -", length(fMAF)))
            myNotes <- c(myNotes, paste("Files successfully processed -", filesFinished))
            myNotes <- c(myNotes, paste("Files dropped due to incomptibility -", numFilesInvalid))

            # send info to log file
            logEvent(eventTitle = myTitle,
                      eventCall = myCall,
                      notes = myNotes)
        }
    }
}

.mafHeaderCheck <- function(file) {
    headValid <- character()
    validMAFheader <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                        "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
                        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
                        "Tumor_Sample_UUID")
    if (! all(read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                     na.strings = " ", nrows=1)[,c(1,5,6,7,8,9,10,11,12,13,16,33)]
        == validMAFheader)) {
        headValid <- c(headValid, "Header of input MAF file did not match the requirements.")
    }
    return(headValid)
}

.rSNVheadercheck <- function(file) {
    headValid <- character()
    validSNVheader <- c("sym", "chr", "start", "end","strand", "class", "type", "aRef",
                        "a1", "a2", "TCGA","UUID")
    if (! all(read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                     na.strings = " ", nrows=1)[,c(1,2,3,4,5,6,7,8,9,10,11,12)]
        == validSNVheader)) {
        headValid <- c(headValid, "Header of input rSNV file did not match the requirements.")
    }
    return(headValid)
}

.fileValidity <- function(file) {
    validExt <- character()
    if (! (strsplit(file, "\\.")[[1]][-1] %in% c("maf", "txt", "tsv"))) {
        validExt <- c(validExt, "Invalid input MAF file extension.")
    }
    return(validExt)
}

# [END]
