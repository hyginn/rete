# importSNV.TCGA.R

#' Import MAF files from source and converts them to rSNV files.
#'
#' \code{importSNV.TCGA} imports list of maf file names, creates rSNV or adds to
#' existing rSNV file, return rSNV file
#'
#' @section MAF file:...
#'
#' @param fMAF vector of MAF file names
#' @param rSNV rSNV file to merged new MAF data with
#' @param na.rm remove all rows with NA values as attributes and drop from DF
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#' @return rSNV file name containg data from all MAF files
#'
#' @family
#'
#'   ## @seealso \code{\link{importX}} importX() is used to import ...
#'
#'   ## @examples ## \dontrun{ ## importSNV.TCGA(IN, OUT) ## }
#' @export
importSNV.TCGA <- function(  fMAF,
                             rSNV = NULL,
                             na.rm = TRUE,
                             silent = FALSE,
                             writeLog = TRUE) {

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
    if (missing(rSNV)) {
        write("# metadata Start\n", "MAFtoSNV", append=FALSE)
        # addMetadata <- function() {}
        write("sym\tchr\tstart\tend\tstrand\tclass\ttype\taRef\ta1\ta2\tTCGA\tUUID", "MAFtoSNV", append=TRUE)
    } else {
        rSNVf <- read.table(rSNV, header = TRUE)[,1]
        rowsBefore <- length(rSNVf)
    }

    myNotes <- c("Number of rows originally in rSNV - ", rowsBefore)

    # ==== READ DATA ===========================================================

    for (i in fMAF) {

        hugo_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE)[,0]
        chr_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE)[,5]
        start_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE)[,6]
        end_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE)[,7]
        strand_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE)[,8]
        var_class_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE)[,9]
        var_type_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE)[,10]
        ref_allele_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE)[,11]
        tum_all1_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE)[,12]
        tum_all2_ <- read.table(i, header = TRUE, fill = NA,  stringsAsFactors = FALSE)[,13]
        tum_samp_barcode_ <- read.table(i, header = TRUE, fill = NA,  stringsAsFactors = FALSE)[,16]
        tum_sam_uuid_ <- read.table(i, header = TRUE, fill = NA,  stringsAsFactors = FALSE)[,33]

        # check if variant_class contains valid values from known list

        classVknown <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                    "Nonsense_Mutation", "Silent", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation",
                    "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR1" , "Intron", "RNA", "Targeted_Region")
        myVclasses <- unique(var_class_)
        if (! all(myVclasses %in% classVknown)) {
            # drop file
            classValid <- FALSE
        } else {
            classValid <- TRUE
        }

        # check if vairant_type contains valid values from known list
        typeVknown <- c("SNP", "DNP", "TNP", "ONP", "INS", "DEL", "Consolidated2")
        myVtypes <- unique(var_type_)
        if (! all(myVtypes %in% typeVknown)) {
            # drop file
            typeValid <- FALSE
        } else {
            typeValid <- TRUE
        }

        # check chromosome value
        typeChr <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
                     "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","X","Y","Z")
        myVchr <- unique(chr_)
        if (! all(myVchr %in% typeChr)) {
            # drop file
            chrValid <- FALSE
        } else {
            chrValid <- TRUE
        }

        # check start position and end are integers
        if (classValid & typeValid & chrValid = FALSE) {
            # there was an invalid value
            # report error in file and log
            # drop file
            # go to next
        } else {
            # all tests passed
            # add data to rSNV
            fileasDT <- cbind(hugo_ , chr_ , start_, end_, strand_, var_class_,
                              var_type_, ref_allele_, tum_all1_, tum_all2_,
                              tum_samp_barcode_, tum_sam_uuid_)
            write.table(fileasDT,  file=rSNV, sep="\t", append=TRUE, row.names =
                            FALSE, col.names = FALSE, quote = FALSE)
        }
    }


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
                 input = character(fMAF, rSNV),
                 notes = myNotes,
                 output = myOutput)
    }

    return(rSNV)
}


# [END]
