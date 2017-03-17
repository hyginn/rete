# importNet.TCGA.R

#' Import MAF files from source and converts them to rSNV files.
#'
#' \code{importSNV.TCGA} imports maf, creates rSNV, returns rSNV
#'
#' @section MAF file:...
#'
#' @param fMAF vector of MAF file names
#' @param rSNV rSNV file to merged new MAF data with
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#' @param na.rm remove all rows with NA values as attributes and drop from DF
#' @return rSNV file name containg data from all MAF files
#'
#' @family
#'
#'   ## @seealso \code{\link{fastMap}} fastMap() is used internally to map ENSP
#'   IDs to gene symbols.
#'
#'   ## @examples ## \dontrun{ ## importNet.STRING(IN, OUT) ## }
#' @export
importSNV.TCGA <- function(  fMAF,
                             rSNV,
                             silent = FALSE,
                             writeLog = TRUE,
                             na.rm = TRUE) {

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

    #  Validate the vector of MAF file names parameter by checking
    #  extension (.maf, .txt, .tsv. gz(?), .zip(?)) and header

    #  Validate rSNV as an existing file
    #  If not, create rSNV file as empty
    #  If rSNV exists, load file
    #  write to myNotes how many lines/rows did rSNV start with

    myNotes <- c("Number of rows originally in rSNV - ")

    # ==== READ DATA ===========================================================

    # reading data will occur in multiple steps to avoid storing too much data on
    # the disk drive and cause overload when handeling large files

    # step 1: go through fMAF input vector list of files and open first file
    read::read_delim()
    # open columns 1,5,6,7,8,9,10,11,12,13,16,33
    # step 2: read lines from file and process data as indicated in "EXTRACT" below
    # typecheck and error handling
    # step 3: save processed data to new/merge to rSNV file
    #         merge temp data to dataframe and delete temp data
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
    # write missing/dropped to log HERE
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
        myCall[2] <- sprintf("fMAF = \"%s\", ", paste(fMAF, collapse=" "))
        myCall[3] <- sprintf("rSNV = \"%s\", ", rSNV)
        myCall[4] <- sprintf("silent = %s, ", as.character(silent))
        myCall[5] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall[6] <- sprintf("na.rm = %s)", as.character(na.rm))
        myCall <- paste0(myCall, collapse = "")

        # indicate input object name(s)
        #   (NA for this importSNV.TCGA())

        # Record progress information
        myNotes <- character()
        myNotes <- c(myNotes, sprintf("Files processed - ", nDone))

        # indicate output object name(s)
        myOutput = c("rSNV")

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
#                input = character(),
                 notes = myNotes,
                 output = myOutput)
    }

    return(rSNV)
}


# [END]
