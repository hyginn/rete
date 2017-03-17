# importNet.TCGA.R

#' Import MAF files from source and coverts them to rSNV files.
#'
#' \code{importSNV.TCGA} imports maf, creates rSNV, returns rSNV
#'
#' @section MAF file:...
#'
#' @param fMAF vector of MAF file names
#' @param rSNV new rSNV file name to be created or the one to be merged with
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param noLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#' @return rSNV file name containg data from all MAF files
#'
#' @family
#'
#'   ## @seealso \code{\link{fastMap}} fastMap() is used internally to map ENSP
#'   IDs to gene symbols.
#'
#'   ## @examples ## \dontrun{ ## importNet.STRING(IN, OUT) ## }
#' @export
importSNV.TCGA <- function(fMAF,
                             rSNV,
                             silent = FALSE,
                             noLog = TRUE) {

    # ==== PARAMETERS ==========================================================

    # NL <- .PlatformLineBreak()
    # cutoffTypes <- c("xS", "xQ", "xN")
    # defaultValues <- c(xS = 950, xQ = 0.9, xN = 10000)
    #
    # if (missing(val)) {
    #     val <- defaultValues[cutoffType]
    # }
    #
    # # Read header and one contents line
    # tmp <- readLines(fName, n = 2)
    #
    # # Parse for requested column. STRING data is " "-delimited.
    # header <- unlist(strsplit(tmp[1], " "))
    # iCol <- which(net == header)  # column that contains the requested net
    #
    # data   <- unlist(strsplit(tmp[2], " "))


    # ==== VALIDATIONS =========================================================

    # General parameter checks
    cR <- character()
    cR <- c(cR, .checkArgs(fMAF,        like = c("FILE_E"),   checkSize = TRUE))
    cR <- c(cR, .checkArgs(rSNV,      like = "FILE_E",        checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,     like = logical(1), checkSize = TRUE))

    if(length(cR) > 0) {
        stop(cR)
    }

    #  Validate the vector of MAF file names parameter
    #  Validate rSNV as an existing file or create new file

    # ==== READ DATA ===========================================================

    # reading data will occur in multiple steps to avoid storing too much data on
    # the disk drive and cause overload when handeling large files

    # step 0: run terminal comand to gzip all folders and subfolders and remove
    #         unwanted files like stdout, stdin, manifest etc..
    # step 1: go through diseases for which data is available and open file
    # step 2: read lines from file and process data as indicated in "EXTRACT" below
    # step 3: save processed data to new/merge to rSNV file
    # step 4: close file
    # step 5: repeat steps 1-4 for next disease MAF file avilable.

    # ==== "EXTRACT" - REQUIRED DATA ========================================================

    # required data from MAF files for rSNV lie in the columns:
    # 1 HUGO SYMBOL
    # 5 chromosome
    # 6 start
    # 7 end
    # 8 strand
    # 9 Variant_Classification
    # 10 Variant_type
    # 11 Reference_allele
    # 12 Tumor_Seq_Allele1
    # 13 Tumor_Seq_Allele2
    # 16 Tumor_Sample_Barcode
    # 33 Tumor_Sample_UUID

    # read each sample from file and parse the columns listed above and write them
    # to the rSNV file.
    # incase of missing value add index of sample entry in rSNV to vector 'missingValues'
    # when we reach the last sample close file and open new file.

    # ==== "HANDLE" - MISSING VALUES =============================================

    # Incase we encounter missing attributes for samples
    # Add the sample entry to the rSNV file but record the line number
    # or the index of that entry.
    # Make a vector of all entries with missing values
    # Report all missing entries in the log Below.

    # can also have and option to drop all entries with missing values

    # ====  =================================================

    # ==== WRITE LOG ===========================================================

    if(writeLog) {

        myTitle <- "importSNV.TCGA"

        # Compile function call record
        myCall <- character()
        myCall[1] <- "importSNV.TCGA("
        myCall[2] <- sprintf("rMAF = \"%s\", ", fName)
        myCall[3] <- sprintf("rSNV = \"%s\", ", net)
        myCall[4] <- sprintf("silent = %s, ", as.character(silent))
        myCall[5] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")

        # indicate input object name(s)
        #   (NA for this importSNV.TCGA())

        # Record progress information
        myNotes <- character()
        myNotes <- c(myNotes, sprintf("Files processed - ", nDone))
        nyNotes <- c(myNotes, sprintf("Files pending - ", nPend))
        myNotes <- c(myNotes, sprintf("Progress in percentage =
                                      Size of Data Processed/Total Size of Data",
                                      nComplete))
        for (i in missingEntry) {
            myNotes <- c(myNotes, sprintf("Found missing entry at index
                                          K of rSNV", K))
        }

        # indicate output object name(s)
        myOutput = c("rSNV")

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
#                input = character(),
                 notes = myNotes,
                 output = myOutput)
    }

    return(gG)
}


# [END]
