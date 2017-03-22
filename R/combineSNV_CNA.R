# combineSNV_CNA.R

#' Combine mutation data from CNA and SNV
#'
#' \code{combineSNV_CNA} loads the vector of preprocessed dataset from filtered CNA or SNV files,
#'  creates hash keys corresponding to the input file type, combine and save the hashed
#'  keys to a gX format file as output file.
#'
#' Details.
#' @section gX files:
#'  gX data frame files contain data of SNV and CNA combined for each gene as one object, sorted per
#'  gene.
#'
#'
#' @param fName A vector of path of fully qualified file names of SNV and/or CNA files.
#' @param fgX The fully qualified filename/path of a gX output file, by default is gX.rds.
#' @param silent Boolean option for writing combiling process information to console, FALSE by default.
#' @param writeLog Boolean option for log results, TRUE by default.
#'
#'
#' @examples
#' \dontrun{combine(fName=file_names, fgX="combined_data.rds", silent=TRUE, writeLog=FALSE)}
#'
#' @export
combineSNV_CNA <- function(fname, fgX = "gX.rds", silent=FALSE, writeLog=TRUE) {
    #' pseudo code
    #'
    #' using hash package, e.g extract keys
    #' if(!require(hash)){
    #'  install.packages("hash")
    #'  library(hash)
    #' }
    #' hashTable <- new.env()
    #' read from input files and store hash
    #' for (i <- 1, i <= length(fname), i++) {
    #'  fcurrent <- fname[i]
    #'  if fcurrent is CNA {
    #'      loadRDS(fcurrent)
    #'      for each variation:
    #'          key = hash(<gene symbol>:CNA:<as.character(round(<copy number>))>)
    #'          if key is in hashTable:
    #'              hashTable[[key]] <- hashTable[[key]] + 1
    #'          else:
    #'              hashTable[[key]] <- 1
    #'  }
    #'  else if fcurrent is MUT {
    #'      loadRDS(fcurrent)
    #'      for each row:
    #'          key = hash(<gene symbol>:SNV:<position>:<variant class>)
    #'          if key is in hashTable:
    #'              hashTable[[key]] <- hashTable[[key]] + 1
    #'          else:
    #'              hashTable[[key]] <- 1
    #' }
    #'
    #'
    #'
    #' retrieve all keys from hashTable
    #' hashVector = hashTable[!is.na(hashTable)] or hashVector <- keys(hashTable)
    #' sort hashVector
    #' strplit them to sym, type, pos and class
    #' store retrieve 4 classes in a dataframe
    #' create NA for not applicable columns in CNAs
    #' for key in hashVector:
    #'  retrieve value from hashTable using keys in hashVector
    #'  insert metadata
    #'
    #' save(dataframe, fgX)
    #'
    #' ToDo add silent and log to read from data files and writing to output files
}

# [END]
