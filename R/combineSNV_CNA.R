# combineSNV_CNA.R

#' Combine mutation data from CNA and SNV
#'
#' \code{combineSNV_CNA} loads the vector of preprossed dataset from filtered CNA or SNV files,
#'  creates hash keys corresponding to the input file type, combine and save the hashed
#'  keys to a gX format file as output file.
#'
#' Details.
#' @section gX files:
#'  gX files contain data of SNV and CNA combined for each gene as one object, sorted per
#'  gene. Tab-delimited text file with header.
#'
#'
#' @param fName A vector of path of fully qualified file names of SNV and/or CNA files.
#' @param fgX The fully qualified filename/path of a gX output file.
#' @param silent Boolean option for writing combiling process information to console.
#' @param noLog Boolean option for log results.
#' @return number of hashed keys stored to the gX output file.
#'
#'
#' @examples
#' \dontrun{combine(fName=file_names, fgX="combined_data.RDS", silent=TRUE, noLog=FALSE)}
#'
#' @export
combineSNV_CNA <- function(fname, fgX, silent=FALSE, noLog=FALSE) {
    #' pseudo code
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
    #'      load(fcurrent)
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
    #' retrieve all keys from hashTable and strplit them to sym, type, pos and class
    #' hashVector = hashTable[!is.na(hashTable)]
    #' sort hashVector
    #' store retrieve 4 classes in a dataframe
    #' for key in hashVector:
    #'  retrieve value from hashTable using keys in hashVector
    #' count <- length(hashVector)
    #'
    #' save(dataframe, fgX)
    #'
    #' ToDo add silent and log to reading from data files and writing to output files
    #' return(count)
}

# [END]
