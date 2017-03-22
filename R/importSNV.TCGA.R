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
                             na.rm = FALSE,
                             silent = FALSE,
                             writeLog = TRUE) {

    myNotes <- character()

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

    meta <- list(type = "rSNV",
                 version = "1.0")

    # ==== Validata fMAF: Extensions  =======================================
    fileValidity <- function(file) {
        if (strsplit(file, "\\.")[[1]][-1] %in% c("maf", "txt", "tsv")) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }

    mafHeaderCheck <- function(file) {
        headValid <- TRUE
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,1] == "Hugo_Symbol") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,5] == "Chromosome") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,6] == "Start_Position") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,7] == "End_Position") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,8] == "Strand") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,9] == "Variant_Classification") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,10] == "Variant_Type") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,11] == "Reference_Allele") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,12] == "Tumor_Seq_Allele1") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,13] == "Tumor_Seq_Allele2") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,16] == "Tumor_Sample_Barcode") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,33] == "Tumor_Sample_UUID") {
            headValid <- FALSE
        }
        return(headValid)
    }

    # ==== Validata rSNV  =======================================

    rSNVheadercheck <- function(file) {
        headValid <- TRUE
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,1] == "sym") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,5] == "chr") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,6] == "start") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,7] == "end") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,8] == "strand") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,9] == "class") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,10] == "type") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,11] == "aRef") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,12] == "a1") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,13] == "a2") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,16] == "TCGA") {
            headValid <- FALSE
        }
        if (! read.table(file, header = FALSE, fill = NA, stringsAsFactors = FALSE,
                         na.strings = " ", nrows=1)[,21] == "UUID") {
            headValid <- FALSE
        }
        return(headValid)
    }

    if (missing(rSNV)) {
        rSNV <- "MAFtoSNV"

        # add metdata
        for (name in names(meta)) {
            attr(rSNV, name) <- meta[[name]]
        }

        write("sym\tchr\tstart\tend\tstrand\tclass\ttype\taRef\ta1\ta2\tTCGA\tUUID",
              rSNV, append=FALSE)

        snvValidity <- TRUE
    } else {
        rSNVf <- read.table(rSNV, header = TRUE)[,1]
        rowsBefore <- length(rSNVf)
        myNotes <- c("Number of rows originally in rSNV - ", rowsBefore)
        snvValidity <- rSNVheadercheck(rSNV)
    }

    # ==== READ DATA ===========================================================
    filesFinished <- 0
    numFilesInvalid <- 0
    for (i in fMAF) {

        if (!silent) {
            .pBar(i, length(fMAF))
        }

        if (fileValidity(i) & mafHeaderCheck(i) & snvValidity) {
            hugo_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE,
                                na.strings = " ")[,1]
            chr_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE,
                               na.strings = " ")[,5]
            start_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE,
                                 na.strings = " ")[,6]
            end_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE,
                               na.strings = " ")[,7]
            strand_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE,
                                  na.strings = " ")[,8]
            var_class_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE,
                                     na.strings = " ")[,9]
            var_type_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE,
                                    na.strings = " ")[,10]
            ref_allele_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE,
                                      na.strings = " ")[,11]
            tum_all1_ <- read.table(i, header = TRUE, fill = NA, stringsAsFactors = FALSE,
                                    na.strings = " ")[,12]
            tum_all2_ <- read.table(i, header = TRUE, fill = NA,  stringsAsFactors = FALSE,
                                    na.strings = " ")[,13]
            tum_samp_barcode_ <- read.table(i, header = TRUE, fill = NA,  stringsAsFactors = FALSE,
                                            na.strings = " ")[,16]
            tum_sam_uuid_ <- read.table(i, header = TRUE, fill = NA,  stringsAsFactors = FALSE,
                                        na.strings = " ")[,33]

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

            # check chromosome value is valid
            typeChr <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
                         "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")
            myVchr <- unique(chr_)
            if (! all(myVchr %in% typeChr)) {
                # drop file
                chrValid <- FALSE
            } else {
                chrValid <- TRUE
            }

            if (classValid & typeValid & chrValid == FALSE) {
                # there was an invalid value
                # report error in file and log
                # drop file
                # go to next
                myNotes <- c(paste(i, " - File was dropped due to invalid attributes in column chr/class/type"))
                numFilesInvalid <- numFilesInvalid + 1
            } else {
                # all tests passed
                # add data to rSNV
                fileasDT <- cbind(hugo_ , chr_ , start_, end_, strand_, var_class_,
                                  var_type_, ref_allele_, tum_all1_, tum_all2_,
                                  tum_samp_barcode_, tum_sam_uuid_)
                if (na.rm) {
                    na.omit(fileasDT)
                }
                write.table(fileasDT,  file=rSNV, sep="\t", append=TRUE, row.names =
                                FALSE, col.names = FALSE, quote = FALSE)
            }
            filesFinished <- filesFinished + 1
        }


        # ==== WRITE LOG ===========================================================

        if(writeLog) {

            myTitle <- "importSNV.TCGA"

            # Compile function call record
            myCall <- character()
            myCall[1] <- "importSNV.TCGA("
            myCall[2] <- sprintf("fMAF = \"%s\", ", paste(fMAF, collapse=" "))
            myCall[3] <- sprintf("rSNV = \"%s\", ", rSNV)
            myCall[3] <- sprintf("silent = %s, ", as.character(silent))
            myCall[4] <- sprintf("writeLog = %s)", as.character(writeLog))
            myCall[5] <- sprintf("na.rm = %s)", as.character(na.rm))
            myCall <- paste0(myCall, collapse = "")

            # indicate input object name(s)
            #   (NA for this importSNV.TCGA())

            # Record progress information
            myNotes <- c(myNotes, sprintf("Files processed - ", filesFinished))
            myNotes <- c(myNotes, sprintf("Total number of files in input fMAF - ", length(fMAF)))
            myNotes <- c(myNotes, sprintf("Files dropped due to incomptibility - "), numFilesInvalid)

            # indicate output object name(s)
            myOutput = c("rSNV")

            # send info to log file
            logEvent(eventTitle = myTitle,
                      eventCall = myCall,
                      #input = character(),
                      notes = myNotes,
                      output = myOutput)
        }
        return(rSNV)
    }
}


# [END]
