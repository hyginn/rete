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
#' \dontrun{combine(fName=file_names, fgX="combined_data.txt", silent=TRUE, noLog=FALSE)}
#'
#' @export
combineSNV_CNA <- function(fname, fgX, silent=FALSE, noLog=FALSE) {
    #' pseudo code
    #' hashVector <- ()
    #' index <- 1
    #' read from input files and store hash
    #' for (i <- 1, i <= length(fname), i++) {
    #'  fcurrent <- fname[i]
    #'  if fcurrent is CNA {
    #'      loadRDS(fcurrent)
    #'      for each variation:
    #'          cnaCount <- count copy number variation
    #'          hashVector[index] <- hash(cnaCount)
    #'          index++
    #'  }
    #'  else if fcurrent is MUT {
    #'      load(fcurrent)
    #'      for each row:
    #'          mutCout <- count observed mutation
    #'          hashVector[index] <- hash(mutCount)
    #'          index++
    #'  }
    #' }
    #' sort hash
    #' hashVector <- sort(hash_vector)
    #' count <- (length(hash_vector))
    #'
    #' write to output gX file from hash
    #' fout <- file(fgX)
    #' for (j <- 1, j <= count, j++) {
    #'  line <- read hash hashVector[j]
    #'  write line to fout
    #' }
    #' close(fout)
    #'
    #' ToDo add silent and log to reading from data files and writing to output files
    #' return(count)
}

# [END]
